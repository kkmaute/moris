/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm.cpp
 *
 */

#include "cl_OPT_Algorithm.hpp"
#include "cl_OPT_Problem.hpp"
#include "HDF5_Tools.hpp"

#include "fn_linspace.hpp"

namespace moris::opt
{

    // -------------------------------------------------------------------------------------------------------------

    Algorithm::Algorithm()
    {
        // set size of FD epsilons to zero to indicate that no values have been set
        mFiniteDifferenceEpsilons.set_size( 0, 0 );

        // set size of ADVS for FD to zero to indicate that no valuves have been set
        mFiniteDifferenceADVs.set_size( 0, 0 );

        // set sensitivity analysis type
        set_sensitivity_analysis_type( mSAType );
    }

    // -------------------------------------------------------------------------------------------------------------

    Algorithm::~Algorithm()
    {
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::compute_design_criteria( Vector< real >& aADVs )
    {
        // check that only processor with rank 0 enters this routine
        MORIS_ERROR( par_rank() == 0,
                "Algorithm::compute_design_criteria - processor with rank > 0 detected.\n" );

        // Coordinate computation of design criteria
        mRunning = opt::Task::compute_criteria_forward_analysis;
        this->communicate_running_status();

        // Compute design criteria
        this->mProblem->compute_design_criteria( aADVs );

        // Set flag that for this design gradients have not been computed yet
        mGradientsHaveBeenComputed = false;
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_objectives()
    {
        // Get objective value(s) from problem
        return mProblem->get_objectives();
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_constraints()
    {
        // Get objective value(s) from problem
        return mProblem->get_constraints();
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::compute_design_criteria_gradients( const Vector< real >& aADVs )
    {
        // Compute gradients of design criteria
        ( this->*compute_design_criteria_gradients_by_type )( aADVs );

        // Set flag that for this design gradients have been computed yet
        mGradientsHaveBeenComputed = true;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::compute_design_criteria_gradients_analytically( const Vector< real >& aADVs )
    {
        // check that only processor with rank 0 enters this routine
        MORIS_ERROR( par_rank() == 0,
                "Algorithm::compute_design_criteria_gradients_analytically - processor with rank > 0 detected.\n" );

        // Coordinate computation of design criteria
        mRunning = opt::Task::compute_criteria_gradients_analytically;
        this->communicate_running_status();

        // Compute design criteria gradients
        mProblem->compute_design_criteria_gradients( aADVs );
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::compute_design_criteria_gradients_fd_fwbw( const Vector< real >& aADVs )
    {
        // Copy ADV vector such that it can be altered
        Vector< real > tADVs = aADVs;

        // number of ADVs with respect to which sensitivities are computed by FD
        uint tNumFDadvs = mFiniteDifferenceADVs.numel();

        // save objective and constraints of unperturbed ADV vector
        // note: need to make a copy as vector in problem will change
        const Matrix< DDRMat > tObjectiveOrg   = this->get_objectives();
        const Matrix< DDRMat > tConstraintsOrg = this->get_constraints();

        // number of objectives
        uint tNumberOfObjectives = tObjectiveOrg.numel();

        // number of constraints
        uint tNumberOfConstraints = tConstraintsOrg.numel();

        // set size of gradient vectors for objective and constraints
        mFDObjectiveGradients.set_size( tNumberOfObjectives, tNumFDadvs );
        mFDObjectiveConstraints.set_size( tNumberOfConstraints, tNumFDadvs );

        // FD each ADV
        for ( uint tIndex = 0; tIndex < tNumFDadvs; tIndex++ )
        {
            // get ADV index
            uint tADVIndex = mFiniteDifferenceADVs( tIndex );

            // Perturb
            tADVs( tADVIndex ) += mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex );

            // Coordinate computation of design criteria
            mRunning = opt::Task::compute_criteria_finite_difference_analysis;
            this->communicate_running_status();

            // Compute design criteria
            this->mProblem->compute_design_criteria( tADVs );

            // get objective and constraints
            const Matrix< DDRMat >& tObjectivePerturbed   = this->get_objectives();
            const Matrix< DDRMat >& tConstraintsPerturbed = this->get_constraints();

            // fw/bw finite difference
            for ( uint tObjIndex = 0; tObjIndex < tNumberOfObjectives; ++tObjIndex )
            {
                mFDObjectiveGradients( tObjIndex, tIndex ) =
                        ( tObjectivePerturbed( tObjIndex ) - tObjectiveOrg( tObjIndex ) ) / mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex );
            }

            for ( uint tConIndex = 0; tConIndex < tNumberOfConstraints; ++tConIndex )
            {
                mFDObjectiveConstraints( tConIndex, tIndex ) =
                        ( tConstraintsPerturbed( tConIndex ) - tConstraintsOrg( tConIndex ) ) / mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex );
            }

            // Restore ADV
            tADVs( tADVIndex ) = aADVs( tADVIndex );
        }

        // Restore objective and gradients
        mProblem->set_objectives_and_constraints( tObjectiveOrg, tConstraintsOrg );
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::compute_design_criteria_gradients_fd_central( const Vector< real >& aADVs )
    {
        // Copy ADV vector such that it can be altered
        Vector< real > tADVs = aADVs;

        // number of ADVs with respect to which sensitivities are computed by FD
        uint tNumFDadvs = mFiniteDifferenceADVs.numel();

        // get a copy of objective and constraints of unperturbed ADV vector
        const Matrix< DDRMat > tObjectiveOrg   = this->get_objectives();
        const Matrix< DDRMat > tConstraintsOrg = this->get_constraints();

        // number of objectives
        uint tNumberOfObjectives = tObjectiveOrg.numel();

        // number of constraints
        uint tNumberOfConstraints = tConstraintsOrg.numel();

        // set size of gradient vectors for objective and constraints
        mFDObjectiveGradients.set_size( tNumberOfObjectives, tNumFDadvs );
        mFDObjectiveConstraints.set_size( tNumberOfConstraints, tNumFDadvs );

        // FD each ADV
        for ( uint tIndex = 0; tIndex < tNumFDadvs; tIndex++ )
        {
            // get ADV index
            uint tADVIndex = mFiniteDifferenceADVs( tIndex );

            // Perturb
            tADVs( tADVIndex ) += mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex );

            // Coordinate computation of design criteria
            mRunning = opt::Task::compute_criteria_finite_difference_analysis;
            this->communicate_running_status();

            // Compute design criteria
            this->mProblem->compute_design_criteria( tADVs );

            // get objective and constraints
            // note: need to make a copy as vector in problem will change
            const Matrix< DDRMat > tObjectivePerturbedPlus   = this->get_objectives();
            const Matrix< DDRMat > tConstraintsPerturbedPlus = this->get_constraints();

            // Perturb
            tADVs( tADVIndex ) -= 2.0 * mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex );

            // Coordinate computation of design criteria
            mRunning = opt::Task::compute_criteria_finite_difference_analysis;
            this->communicate_running_status();

            // Compute design criteria
            this->mProblem->compute_design_criteria( tADVs );

            // get objective and constraints; note: do not need to make a copy
            const Matrix< DDRMat >& tObjectivePerturbedMinus   = this->get_objectives();
            const Matrix< DDRMat >& tConstraintsPerturbedMinus = this->get_constraints();

            // central finite difference
            for ( uint tObjIndex = 0; tObjIndex < tNumberOfObjectives; ++tObjIndex )
            {
                mFDObjectiveGradients( tObjIndex, tIndex ) =
                        ( tObjectivePerturbedPlus( tObjIndex ) - tObjectivePerturbedMinus( tObjIndex ) ) / ( 2.0 * mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex ) );
            }

            for ( uint tConIndex = 0; tConIndex < tNumberOfConstraints; ++tConIndex )
            {
                mFDObjectiveConstraints( tConIndex, tIndex ) =
                        ( tConstraintsPerturbedPlus( tConIndex ) - tConstraintsPerturbedMinus( tConIndex ) ) / ( 2.0 * mFiniteDifferenceEpsilons( tADVIndex, mPerturbationSizeIndex ) );
            }

            // Restore ADV
            tADVs( tADVIndex ) = aADVs( tADVIndex );
        }

        // Restore objective and gradients
        mProblem->set_objectives_and_constraints( tObjectiveOrg, tConstraintsOrg );
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_objective_gradients()
    {
        // Check that gradients have been computed
        MORIS_ERROR( mGradientsHaveBeenComputed,
                "Algorithm::get_objective_gradients - objective gradients not available for current advs.\n" );

        // return gradients
        return ( this->*get_objective_gradients_by_type )();
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_objective_gradients_analytically()
    {
        return mProblem->get_objective_gradients();
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_objective_gradients_by_fd()
    {
        return mFDObjectiveGradients;
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_constraint_gradients()
    {
        // Check that gradients have been computed
        MORIS_ERROR( mGradientsHaveBeenComputed,
                "AAlgorithm::get_constraint_gradients - constraint gradients not available for current advs.\n" );

        // return gradients
        return ( this->*get_constraint_gradients_by_type )();
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_constraint_gradients_analytically()
    {
        return mProblem->get_constraint_gradients();
    }

    //----------------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Algorithm::get_constraint_gradients_by_fd()
    {
        return mFDObjectiveConstraints;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::set_sensitivity_analysis_type( SA_Type aType )
    {
        // Set gradient function pointers
        switch ( aType )
        {
            case opt::SA_Type::analytical:
            {
                compute_design_criteria_gradients_by_type = &Algorithm::compute_design_criteria_gradients_analytically;
                get_objective_gradients_by_type           = &Algorithm::get_objective_gradients_analytically;
                get_constraint_gradients_by_type          = &Algorithm::get_constraint_gradients_analytically;
                return;
            }
            case opt::SA_Type::backward:
            {
                mFiniteDifferenceEpsilons = -1.0 * mFiniteDifferenceEpsilons;
                [[fallthrough]];
            }
            case opt::SA_Type::forward:
            {
                compute_design_criteria_gradients_by_type = &Algorithm::compute_design_criteria_gradients_fd_fwbw;
                get_objective_gradients_by_type           = &Algorithm::get_objective_gradients_by_fd;
                get_constraint_gradients_by_type          = &Algorithm::get_constraint_gradients_by_fd;
                break;
            }
            case opt::SA_Type::central:
            {
                compute_design_criteria_gradients_by_type = &Algorithm::compute_design_criteria_gradients_fd_central;
                get_objective_gradients_by_type           = &Algorithm::get_objective_gradients_by_fd;
                get_constraint_gradients_by_type          = &Algorithm::get_constraint_gradients_by_fd;
                break;
            }
            default:
            {
                MORIS_ERROR( false, "Algorithm::set_finite_differencing - Sensitivity analysis type not implemented.\n" );
            }
        }

        // this point is only reached when FD is used; thus FD schemes are initialized
        this->initialize_finite_difference_schemes();
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::set_finite_difference_perturbation_size_index( uint aPerturbationSizeIndex )
    {
        mPerturbationSizeIndex = aPerturbationSizeIndex;
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Algorithm::set_finite_difference_advs( const Matrix< DDUMat >& aFiniteDifferenceADVs )
    {
        mFiniteDifferenceADVs = aFiniteDifferenceADVs;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::communicate_running_status()
    {
        // check that incoming status of processors with rank larger 1 is "wait"
        MORIS_ERROR( par_rank() > 0 ? mRunning == Task::wait : true,
                "Algorithm::communicate_running_status - incoming status of processors with rank larger 1 should be wait\n" );

        // convert enum into uint for MPI communication
        uint tMsg = (uint)mRunning;

        // processor 0 sends running status to all other processors
        broadcast( tMsg );

        // update running status
        mRunning = (Task)tMsg;

        // Increase optimization iteration counter if forward analysis
        if ( mRunning == Task::compute_criteria_forward_analysis )
        {
            MORIS_LOG_ITERATION();
        }
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm::dummy_solve()
    {
        // Communicate that these procs need to start running
        this->communicate_running_status();

        // Create dummy ADVs
        Vector< real > tDummyADVs;

        // Perform requested analysis type until exit status received
        while ( mRunning != Task::exit )
        {
            switch ( mRunning )
            {
                case Task::exit:
                {
                    // Do nothing on exit
                    break;
                }
                case Task::wait:
                {
                    // Do nothing on exit
                    break;
                }
                case Task::compute_criteria_forward_analysis:
                case Task::compute_criteria_finite_difference_analysis:
                {
                    // Compute design criteria
                    this->mProblem->compute_design_criteria( tDummyADVs );
                    break;
                }
                case Task::compute_criteria_gradients_analytically:
                {
                    // Compute analytically derivatives of design criteria
                    this->mProblem->compute_design_criteria_gradients( tDummyADVs );
                    break;
                }
                default:
                {
                    MORIS_ERROR( false,
                            "Algorithm::dummy_solve - undefined running status.\n" );
                }
            }

            // set task to wait
            mRunning = Task::wait;

            // Communicate running status so these processors know when to exit
            this->communicate_running_status();
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Algorithm::write_advs_to_file( const Vector< real >& aADVs )
    {
        // Get iteration from global clock
        uint tOptIter = gLogger.get_opt_iteration();

        // Create file name
        std::string tRestartFileName =
                "ADV_Alg_" + std::to_string( mCurrentOptAlgInd ) + "_Iter_" + std::to_string( tOptIter ) + ".hdf5";

        // Create file
        // Note: only processor 0 creates this file; therefore no parallel name extension is used
        hid_t tFileID = create_hdf5_file( tRestartFileName, false );

        // Write advs and upper/lower bounds to file
        herr_t tStatus = 0;

        save_vector_to_hdf5_file( tFileID, "ADVs", aADVs.data(), tStatus );
        save_vector_to_hdf5_file( tFileID, "UpperBounds", mProblem->get_upper_bounds().data(), tStatus );
        save_vector_to_hdf5_file( tFileID, "LowerBounds", mProblem->get_lower_bounds().data(), tStatus );

        if ( mProblem->get_ijklIDs().numel() == aADVs.size() )
        {
            save_matrix_to_hdf5_file( tFileID, "IjklIds", mProblem->get_ijklIDs(), tStatus );
        }

        save_scalar_to_hdf5_file( tFileID, "NumProcs", par_size(), tStatus );

        // Close file
        close_hdf5_file( tFileID );
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Algorithm::initialize_finite_difference_schemes()
    {
        // get number of advs
        uint tNumberOfAdvs = mProblem->get_num_advs();

        // determine with respect to which advs sensitivities are compute by FD; if list not set in input file
        // all sensitivities with respect to all advs are computed
        uint tNumFDadvs = mFiniteDifferenceADVs.numel();

        if ( tNumFDadvs == 0 )
        {
            // Set number of FD-ADVs to number of ADVs
            tNumFDadvs = tNumberOfAdvs;

            // fill adv list with all advs
            mFiniteDifferenceADVs = linspace< uint >( 0, tNumberOfAdvs - 1, tNumberOfAdvs );
        }
        else
        {
            // check that requested adv indices exists
            for ( uint tIndex = 0; tIndex < tNumFDadvs; ++tIndex )
            {
                MORIS_ERROR( mFiniteDifferenceADVs( tIndex ) < tNumberOfAdvs,
                        "Algorithm::initialize_finite_difference_schemes - %s",
                        "Requested Adv index for FD too large.\n" );
            }
        }

        // Finite differencing perturbation size
        if ( mFiniteDifferenceEpsilons.numel() != 0 )
        {
            // check if only one vector of perturbation sizes is provided
            if ( mFiniteDifferenceEpsilons.n_rows() == 1 )
            {
                uint tNumEvals = mFiniteDifferenceEpsilons.n_cols();
                mFiniteDifferenceEpsilons.resize( tNumberOfAdvs, tNumEvals );

                // assign same perturbation values to all ADVs
                for ( uint tIndex = 1; tIndex < tNumberOfAdvs; tIndex++ )
                {
                    mFiniteDifferenceEpsilons( { tIndex, tIndex }, { 0, tNumEvals - 1 } ) =
                            mFiniteDifferenceEpsilons( { 0, 0 }, { 0, tNumEvals - 1 } );
                }
            }

            // check that matrix of perturbation sizes has correct size
            MORIS_ERROR( mFiniteDifferenceEpsilons.n_rows() == tNumberOfAdvs,
                    "Algorithm::initialize_finite_difference_schemes - %s",
                    "Number of rows in finite_difference_epsilons must match the number of ADVs." );
        }
        else
        {
            MORIS_ERROR( false,
                    "Algorithm::initialize_finite_difference_schemes - %s",
                    "At least one value for the FD perturbation size needs to be set.\n" );
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Algorithm::set_restart_index( uint aRestartIndex )
    {
        mRestartIndex = aRestartIndex;

        // calculating updated max its
        mMaxIterations = mMaxIterationsInitial - mRestartIndex;
    }

    // -------------------------------------------------------------------------------------------------------------

}    // namespace moris::opt
