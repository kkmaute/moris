/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem.cpp
 *
 */

#include "cl_OPT_Problem.hpp"
#include "op_plus.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"
#include "HDF5_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "fn_stringify_matrix.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------
        // Public functions
        // -------------------------------------------------------------------------------------------------------------

        Problem::Problem(
                ParameterList&                         aParameterList,
                std::shared_ptr< Criteria_Interface >& aInterface )
        {
            // Set interface
            mInterface = aInterface;

            // Parameters: restart file name
            mRestartFile = aParameterList.get< std::string >( "restart_file" );
        }

        // -------------------------------------------------------------------------------------------------------------

        Problem::~Problem()
        {
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::initialize()
        {
            // Trace optimization problem
            Tracer tTracer( "OPT", "Problem", "Initialize" );

            // Initialize ADVs
            mInterface->initialize( mADVs, mLowerBounds, mUpperBounds, mIjklIds );

            // Override interface ADVs
            this->override_advs();

            // Log number of optimization variables
            MORIS_LOG_SPEC( "Number of optimization variables", mADVs.numel() );

            // Check for proper dimensions of ADV vector and its upper and lower bound vectors
            MORIS_ERROR( mADVs.numel() > 0 ? mADVs.n_cols() == 1 : true,
                    "ADVs vector needs to be column vector.\n" );
            MORIS_ERROR( mLowerBounds.numel() > 0 ? mLowerBounds.n_cols() == 1 : true,
                    "ADV lower bound vector needs to be column vector.\n" );
            MORIS_ERROR( mUpperBounds.numel() > 0 ? mUpperBounds.n_cols() == 1 : true,
                    "ADV lower bound vector needs to be column vector.\n" );

            MORIS_ERROR( mADVs.n_rows() == mLowerBounds.n_rows() and mADVs.n_rows() == mUpperBounds.n_rows(),
                    "ADVs and its lower and upper bound vectors need to have same length.\n" );

            // Read ADVs from restart file
            if ( par_rank() == 0 && !mRestartFile.empty() )
            {
                this->read_restart_file();
            }

            // Initialize constraint types
            mConstraintTypes = this->get_constraint_types();

            MORIS_ERROR( this->check_constraint_order(), "The constraints are not ordered properly (Eq, Ineq)" );    // TODO move call to alg

            // Set number of objectives and constraints
            mNumObjectives  = 1;
            mNumConstraints = mConstraintTypes.numel();

            // Log number of objectives and constraints
            MORIS_LOG_SPEC( "Number of objectives ", mNumObjectives );
            MORIS_LOG_SPEC( "Number of constraints", mNumConstraints );
        }

        // -------------------------------------------------------------------------------------------------------------

        bool
        Problem::check_constraint_order()
        {
            // Checks that the constraint order is  1. equality constraints   (typ=0)
            //                                      2. inequality constraints (typ=1)

            int tSwitches = 1;
            for ( uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++ )
            {
                if ( mConstraintTypes( tConstraintIndex ) == tSwitches )
                {
                    tSwitches -= 1;
                }
            }
            return ( tSwitches >= 0 );
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Problem::get_objectives()
        {
            MORIS_ERROR( par_rank() == 0,
                    "Problem::get_objectives - called by another processor but zero.\n" );

            // log objective
            MORIS_LOG_SPEC( "Objective", ios::stringify_log( mObjectives ) );

            return mObjectives;
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Problem::get_constraints()
        {
            MORIS_ERROR( par_rank() == 0,
                    "Problem::get_constraints - called by another processor but zero.\n" );

            // log constraints
            MORIS_LOG_SPEC( "Constraints", ios::stringify_log( mConstraints ) );

            return mConstraints;
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::set_objectives_and_constraints(
                const Matrix< DDRMat >& aObjectives,
                const Matrix< DDRMat >& aConstraints )
        {
            // check for proper size of input matrices
            MORIS_ERROR( aObjectives.n_rows() == mObjectives.n_rows() and aObjectives.n_rows() == mObjectives.n_rows(),
                    "Problem::set_objectives_and_constraints - input objective vector has incorrect size\n." );

            MORIS_ERROR( aConstraints.n_rows() == mConstraints.n_rows() and aConstraints.n_rows() == mConstraints.n_rows(),
                    "Problem::set_objectives_and_constraints - input objective vector has incorrect size\n." );

            // copy objectives and constraints
            mObjectives  = aObjectives;
            mConstraints = aConstraints;
        }

        // -------------------------------------------------------------------------------------------------------------

        uint
        Problem::get_num_equality_constraints()
        {
            uint tNumEqualityConstraints = 0;

            // Loop over the constraint types matrix
            for ( uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++ )
            {
                if ( mConstraintTypes( tConstraintIndex ) == 0 )
                {
                    tNumEqualityConstraints++;
                }
            }

            return tNumEqualityConstraints;
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::compute_design_criteria( Matrix< DDRMat >& aNewADVs )
        {
            // check for proper size of input ADV vector (on processor 0 only as all other processors
            // have zero-size ADV vector
            if ( par_rank() == 0 )
            {
                // Check that input ADV vector has correct size
                MORIS_ERROR( mADVs.n_rows() == aNewADVs.n_rows() && mADVs.n_cols() == aNewADVs.n_cols(),
                        "Problem::compute_design_criteria - size of ADV vectors does not match.\n" );
            }

            // compute design criteria
            mCriteria = mInterface->get_criteria( aNewADVs );

            // compute objective and constraints (on processor 0 only)
            if ( par_rank() == 0 )
            {
                // save the new ADV vector ( this is done after the criteria is computed since get-criteria might make some changes
                //  in the adv vector e.g reinitialization of the ADV vector)
                mADVs = aNewADVs;

                // log criteria and ADVs
                MORIS_LOG_SPEC( "Criteria", ios::stringify_log( mCriteria ) );

                if ( mADVs.numel() > 0 )
                {
                    MORIS_LOG_SPEC( "MinADV", mADVs.min() );
                    MORIS_LOG_SPEC( "MaxADV", mADVs.max() );
                }

                // compute objective
                mObjectives = this->compute_objectives();

                MORIS_ASSERT( mObjectives.numel() == mNumObjectives,
                        "compute_design_criteria - Number of objectives is incorrect (%zu vs %i).\n",
                        mObjectives.numel(),
                        mNumObjectives );

                // compute constraints
                mConstraints = this->compute_constraints();

                MORIS_ASSERT( mConstraints.numel() == mNumConstraints,
                        "compute_design_criteria - Number of constraints is incorrect (%zu vs %i).\n",
                        mConstraints.numel(),
                        mNumConstraints );
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::compute_design_criteria_gradients( const Matrix< DDRMat >& aNewADVs )
        {
            // Check that input ADV vector is identical to the one stored in problem
            // as otherwise new criteria evaluation would be needed; this is done only
            // on processor 0 only all other processors as have zero-size ADV vector
            if ( par_rank() == 0 )
            {
                MORIS_ERROR( norm( mADVs - aNewADVs ) <= mADVNormTolerance * norm( mADVs ),
                        "Problem::compute_design_criteria_gradients - given ADV vector does not match ADVs used to compute design criteria.\n" );
            }

            // compute derivatives of design criteria
            mInterface->get_dcriteria_dadv();

            // compute objective and constraints (on processor 0 only)
            if ( par_rank() == 0 )
            {
                // compute gradients of objective
                mObjectiveGradient =
                        this->compute_dobjective_dadv() +    //
                        this->compute_dobjective_dcriteria() * mInterface->get_dcriteria_dadv();

                MORIS_ASSERT( mObjectiveGradient.n_rows() == mNumObjectives and mObjectiveGradient.n_cols() == mADVs.numel(),
                        "Problem::compute_design_criteria_gradients  - gradient of objective matrix has incorrect size.\n." );

                // compute gradients of constraints
                mConstraintGradient =
                        this->compute_dconstraint_dadv() +    //
                        this->compute_dconstraint_dcriteria() * mInterface->get_dcriteria_dadv();

                MORIS_ASSERT( mConstraintGradient.n_rows() == mNumConstraints and mConstraintGradient.n_cols() == mADVs.numel(),
                        "Problem::compute_design_criteria_gradients  - gradient of constraint matrix has incorrect size.\n." );
            }
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Problem::get_objective_gradients()
        {
            MORIS_ERROR( par_rank() == 0,
                    "Problem::get_objective_gradients - called by another processor but zero.\n" );

            return mObjectiveGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        const Matrix< DDRMat >&
        Problem::get_constraint_gradients()
        {
            MORIS_ERROR( par_rank() == 0,
                    "Problem::get_constraint_gradients - called by another processor but zero.\n" );

            return mConstraintGradient;
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::scale_solution()
        {
            // TODO Need to decide on a framework for the scaling of the solution
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::update_problem()
        {
            // TODO Need to decide on a framework for the update of the problem
        }

        // -------------------------------------------------------------------------------------------------------------

        void
        Problem::read_restart_file()
        {
            MORIS_LOG_INFO( "Reading ADVs from restart file: %s", mRestartFile.c_str() );

            // Open open restart file
            // Note: only processor 0 reads this file; therefore no parallel file name extension is used
            hid_t tFileID = open_hdf5_file( mRestartFile, false );

            // Define matrix in which to read restart ADVs
            Matrix< DDRMat > tRestartADVs;
            Matrix< DDRMat > tRestartUpperBounds;
            Matrix< DDRMat > tRestartLowerBounds;
            Matrix< IdMat >  tIjklIdsFile;
            sint             tParSizeFile = 0;

            // Read ADVS from restart file
            herr_t tStatus = 0;
            load_matrix_from_hdf5_file( tFileID, "ADVs", tRestartADVs, tStatus );
            load_matrix_from_hdf5_file( tFileID, "UpperBounds", tRestartUpperBounds, tStatus );
            load_matrix_from_hdf5_file( tFileID, "LowerBounds", tRestartLowerBounds, tStatus );

            load_scalar_from_hdf5_file( tFileID, "NumProcs", tParSizeFile, tStatus );

            // FIXME this function also should work if you want to restart with same number procs but new proc layout
            // in this case this if criteria should be deleted
            if ( par_size() != tParSizeFile )
            {
                load_matrix_from_hdf5_file( tFileID, "IjklIds", tIjklIdsFile, tStatus );

                MORIS_ERROR( tIjklIdsFile.numel() == tRestartADVs.numel(), "restart adv vector and ijkl ID vector not of same size" );
                MORIS_ERROR( tIjklIdsFile.numel() == mIjklIds.numel(), "restart ijkl ID vector and new ijkl ID vector not of same size" );
            }

            // Close restart file
            close_hdf5_file( tFileID );

            // Check for matching sizes
            MORIS_ERROR( tRestartADVs.numel() == mADVs.numel(),
                    "Number of restart ADVS does not match.\n" );
            MORIS_ERROR( tRestartUpperBounds.numel() == mUpperBounds.numel(),
                    "Number of restart upper bounds of ADVs does not match.\n" );
            MORIS_ERROR( tRestartLowerBounds.numel() == mLowerBounds.numel(),
                    "Number of restart lower bounds of ADVs does not match.\n" );

            MORIS_LOG_INFO( "Norm of ADV vector - before loading restart %e   after %e.\n",
                    norm( mADVs ),
                    norm( tRestartADVs ) );

            // reorder adv for new proc layout. FIXME see above
            if ( par_size() != tParSizeFile )
            {
                moris::map< moris_id, sint > tIjklIdToPosMap;
                for ( uint Ik = 0; Ik < mIjklIds.numel(); Ik++ )
                {
                    if ( mIjklIds( Ik ) != gNoID )
                    {
                        tIjklIdToPosMap[ mIjklIds( Ik ) ] = Ik;
                    }
                }

                Matrix< DDRMat > tRestartADVsTemp;
                Matrix< DDRMat > tRestartUpperBoundsTemp;
                Matrix< DDRMat > tRestartLowerBoundsTemp;

                tRestartADVsTemp        = tRestartADVs;
                tRestartUpperBoundsTemp = tRestartUpperBounds;
                tRestartLowerBoundsTemp = tRestartLowerBounds;

                tRestartADVs.fill( MORIS_REAL_MAX );
                tRestartUpperBounds.fill( MORIS_REAL_MAX );
                tRestartLowerBounds.fill( MORIS_REAL_MAX );

                for ( uint Ik = 0; Ik < tIjklIdsFile.numel(); Ik++ )
                {
                    if ( tIjklIdToPosMap.key_exists( tIjklIdsFile( Ik ) ) )
                    {
                        sint tIndxed                   = tIjklIdToPosMap.find( tIjklIdsFile( Ik ) );
                        tRestartADVs( tIndxed )        = tRestartADVsTemp( Ik );
                        tRestartUpperBounds( tIndxed ) = tRestartUpperBoundsTemp( Ik );
                        tRestartLowerBounds( tIndxed ) = tRestartLowerBoundsTemp( Ik );
                    }
                    else
                    {
                        tRestartADVs( Ik )        = tRestartADVsTemp( Ik );
                        tRestartUpperBounds( Ik ) = tRestartUpperBoundsTemp( Ik );
                        tRestartLowerBounds( Ik ) = tRestartLowerBoundsTemp( Ik );
                    }
                }

                MORIS_ERROR( tRestartADVs.max() != MORIS_REAL_MAX, "Restart ADV is MORIS_REAL_MAX" );
                MORIS_ERROR( tRestartUpperBounds.max() != MORIS_REAL_MAX, "Restart upper bound is MORIS_REAL_MAX" );
                MORIS_ERROR( tRestartLowerBounds.max() != MORIS_REAL_MAX, "Restart lower bound is MORIS_REAL_MAX" );
            }

            // Copy restart vectors on member variables
            mADVs        = tRestartADVs;
            mUpperBounds = tRestartUpperBounds;
            mLowerBounds = tRestartLowerBounds;

            mRestartFile = "";
        }
    }    // namespace opt
}    // namespace moris
