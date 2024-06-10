/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_Sweep.cpp
 *
 */

#include "cl_OPT_Algorithm_Sweep.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_sum.hpp"
#include "HDF5_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris
{
    namespace opt
    {

        // -------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::Algorithm_Sweep( Parameter_List aParameterList )
        {
            // define which quantities are evaluated
            mEvaluateObjectives  = aParameterList.get< bool >( "evaluate_objectives" );
            mEvaluateConstraints = aParameterList.get< bool >( "evaluate_constraints" );

            mEvaluateObjectiveGradients  = aParameterList.get< bool >( "evaluate_objective_gradients" );
            mEvaluateConstraintGradients = aParameterList.get< bool >( "evaluate_constraint_gradients" );

            // define sampling strategy
            mIncludeBounds    = aParameterList.get< bool >( "include_bounds" );
            mNumEvaluations   = string_to_mat< DDUMat >( aParameterList.get< std::string >( "num_evaluations_per_adv" ) );
            mEvaluationPoints = string_to_mat< DDRMat >( aParameterList.get< std::string >( "custom_adv_evaluations" ) );

            // define finite difference strategy
            mFiniteDifferenceType     = aParameterList.get< std::string >( "finite_difference_type" );
            mFiniteDifferenceEpsilons = string_to_mat< DDRMat >( aParameterList.get< std::string >( "finite_difference_epsilons" ) );
            mFiniteDifferenceADVs     = string_to_mat< DDUMat >( aParameterList.get< std::string >( "finite_difference_adv_indices" ) );

            // define output options
            mSave  = aParameterList.get< bool >( "save" );
            mPrint = aParameterList.get< bool >( "print" );

            // open HDF5 file
            if ( par_rank() == 0 )
            {
                mFileID = create_hdf5_file( aParameterList.get< std::string >( "hdf5_path" ) );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_Sweep::~Algorithm_Sweep()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Algorithm_Sweep::solve(
                uint                       aCurrentOptAlgInd,
                std::shared_ptr< Problem > aOptProb )
        {
            // Trace optimization
            Tracer tTracer( "OptimizationAlgorithm", "Sweep", "Solve" );

            // Initialize
            mCurrentOptAlgInd = aCurrentOptAlgInd;    // set index of current optimization algorithm
            mProblem          = aOptProb;             // set the member variable mProblem to aOptProb

            // Solve optimization problem
            if ( par_rank() == 0 )
            {
                // Run sweep algorithm
                this->sweep_solve();

                // Communicate that optimization has finished
                mRunning = opt::Task::exit;

                this->communicate_running_status();
            }
            else
            {
                // Don't print from these processors
                mPrint = false;

                // Run dummy solve
                this->dummy_solve();
            }

            // update aOptProb
            aOptProb = mProblem;

            return 0;
        }

        //----------------------------------------------------------------------------------------------------------

        void
        Algorithm_Sweep::set_up_evaluation_points()
        {
            // Get number of ADVs
            uint tNumADVs = mProblem->get_num_advs();

            uint tTotalEvaluations = mEvaluationPoints.n_cols();

            if ( mEvaluationPoints.numel() == 0 )
            {
                // Check user input
                if ( mNumEvaluations.numel() == 1 )    // check for global number of evaluations
                {
                    mNumEvaluations.set_size( tNumADVs, 1, mNumEvaluations( 0 ) );

                    // check that number of evaluation points is less than MORIS_UINT_MAX
                    MORIS_ERROR( std::pow( mNumEvaluations( 0 ), tNumADVs ) < (real)MORIS_UINT_MAX,
                            "Algorithm_Sweep::set_up_evaluation_points - Number of evaluations exceeds MORIS_UINT_MAX." );
                }
                else    // check for number of evaluations given per ADV
                {
                    MORIS_ERROR( mNumEvaluations.numel() == tNumADVs,
                            "Must give single number of evaluations for all ADVs or one per ADV, or provide custom evaluation points" );
                }

                // Check lower and upper bounds for equality
                Vector< real > tADVs        = mProblem->get_advs();
                Vector< real > tLowerBounds = mProblem->get_lower_bounds();
                Vector< real > tUpperBounds = mProblem->get_upper_bounds();

                for ( uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++ )
                {
                    // If 1 evaluation, set lower and upper bounds to be equal
                    if ( mNumEvaluations( tADVIndex ) == 1 )
                    {
                        tLowerBounds( tADVIndex ) = tADVs( tADVIndex );
                        tUpperBounds( tADVIndex ) = tADVs( tADVIndex );
                    }

                    // If lower and upper bounds are equal, set 1 evaluation
                    if ( std::abs( tLowerBounds( tADVIndex ) - tUpperBounds( tADVIndex ) ) < MORIS_REAL_EPS )
                    {
                        mNumEvaluations( tADVIndex ) = 1;
                    }
                }

                // Change lower bounds for sweep based on parameter and initialize ADvs
                if ( !mIncludeBounds )
                {
                    for ( uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++ )
                    {
                        tLowerBounds( tADVIndex ) += ( tUpperBounds( tADVIndex ) - tLowerBounds( tADVIndex ) ) / ( mNumEvaluations( tADVIndex ) + 1 );
                    }
                }

                // Set up evaluations
                tTotalEvaluations = 1;
                for ( uint ind = 0; ind < mNumEvaluations.numel(); ind++ )
                {
                    tTotalEvaluations *= mNumEvaluations( ind );
                }

                MORIS_ERROR( tTotalEvaluations > 0,
                        "Algorithm_Sweep::set_up_evaluation_points - Number of evaluations points need to be larger than zero." );

                MORIS_LOG_INFO( "Number of evaluation points: %d", tTotalEvaluations );

                // allocate matrix of evaluation points
                mEvaluationPoints.set_size( tNumADVs, tTotalEvaluations );

                tADVs = tLowerBounds;
                Matrix< DDUMat > tCurrentEvaluations( tNumADVs, 1, 0 );

                // Construct evaluation points
                for ( uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++ )
                {
                    if ( tNumADVs > 0 )
                    {
                        // Assign ADVs
                        for ( uint tADVIndex = 0; tADVIndex < tNumADVs; tADVIndex++ )
                        {
                            mEvaluationPoints( tADVIndex, tEvaluationIndex ) = tADVs( tADVIndex );
                        }

                        // Update ADVs
                        tADVs( 0 ) += ( tUpperBounds( 0 ) - tLowerBounds( 0 ) ) / ( mNumEvaluations( 0 ) + 1 - ( 2 * mIncludeBounds ) );

                        tCurrentEvaluations( 0 ) += 1;

                        for ( uint tADVIndex = 0; tADVIndex < tNumADVs - 1; tADVIndex++ )
                        {
                            if ( tCurrentEvaluations( tADVIndex ) == mNumEvaluations( tADVIndex ) )
                            {
                                // Reset this ADV to the lower bound and increment next ADV
                                tADVs( tADVIndex )               = tLowerBounds( tADVIndex );
                                tCurrentEvaluations( tADVIndex ) = 0;

                                tADVs( tADVIndex + 1 ) +=    //
                                        ( tUpperBounds( tADVIndex ) - tLowerBounds( tADVIndex ) ) / ( mNumEvaluations( tADVIndex ) + 1 - ( 2 * mIncludeBounds ) );
                                tCurrentEvaluations( tADVIndex + 1 ) += 1;
                            }
                        }
                    }
                }
            }
            else
            {
                MORIS_ERROR( mEvaluationPoints.n_rows() == tNumADVs,
                        "Number of rows in custom_adv_evaluations must match the number of ADVs (%d).", tNumADVs );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_Sweep::sweep_solve()
        {
            // Set up evaluation points
            this->set_up_evaluation_points();

            // Initialize FD schemes
            this->initialize_finite_difference_schemes();

            // Total number of perturbation sizes
            uint tTotalEpsilons = mFiniteDifferenceEpsilons.n_cols();

            // Open file and write ADVs/epsilons
            herr_t tStatus = 0;
            if ( mSave )
            {
                moris::save_matrix_to_hdf5_file( mFileID, "adv_evaluations", mEvaluationPoints, tStatus );
                moris::save_matrix_to_hdf5_file( mFileID, "epsilons", mFiniteDifferenceEpsilons, tStatus );
            }

            // Print ADVs/epsilons to be evaluated
            if ( mPrint )
            {
                moris::print( mEvaluationPoints, "adv_evaluations" );
                moris::print( mFiniteDifferenceEpsilons, "epsilons" );
            }

            // Evaluation string
            std::string tEvaluationName;

            // Number of evaluations
            uint tTotalEvaluations = mEvaluationPoints.n_cols();

            // Loop through evaluation points
            for ( uint tEvaluationIndex = 0; tEvaluationIndex < tTotalEvaluations; tEvaluationIndex++ )
            {
                // get the evaluation point
                Vector< real > tEvaluationPoint( mEvaluationPoints.n_rows() );
                for ( uint iADVIndex = 0; iADVIndex < tEvaluationPoint.size(); iADVIndex++ )
                {
                    tEvaluationPoint( iADVIndex ) = mEvaluationPoints( iADVIndex, tEvaluationIndex );
                }

                // Compute design criteria at current evaluation point
                this->compute_design_criteria( tEvaluationPoint );

                // Set evaluation name
                tEvaluationName = " eval_" + std::to_string( tEvaluationIndex + 1 ) + "-" + std::to_string( tTotalEvaluations );

                // Evaluate objectives and constraints
                this->output_objectives_constraints( tEvaluationName );

                // Compute gradients of design criteria at current evaluation point
                if ( mEvaluateObjectiveGradients or mEvaluateConstraintGradients )
                {
                    // Get analytical gradients if requested
                    if ( mFiniteDifferenceType == "none" or mFiniteDifferenceType == "all" )
                    {
                        this->set_sensitivity_analysis_type( SA_Type::analytical );
                        this->compute_design_criteria_gradients( tEvaluationPoint );
                        this->evaluate_objective_gradients( tEvaluationName + " analytical" );
                        this->evaluate_constraint_gradients( tEvaluationName + " analytical" );
                    }

                    // Perform finite difference sensitivity analysis for different perturbation sizes
                    for ( uint tEpsilonIndex = 0; tEpsilonIndex < tTotalEpsilons; tEpsilonIndex++ )
                    {
                        // Reset evaluation name with epsilon data
                        tEvaluationName = " eval_" + std::to_string( tEvaluationIndex + 1 ) + "-" + std::to_string( tTotalEvaluations );
                        if ( mFiniteDifferenceType != "none" )
                        {
                            tEvaluationName += " epsilon_" + std::to_string( tEpsilonIndex + 1 ) + "-" + std::to_string( tTotalEpsilons );
                        }

                        // Compute and/or save gradients based on finite differencing requested
                        if ( mFiniteDifferenceType == "all" )
                        {
                            // set perturbation size index
                            this->set_finite_difference_perturbation_size_index( tEpsilonIndex );

                            // Forward
                            this->set_sensitivity_analysis_type( SA_Type::forward );
                            this->compute_design_criteria_gradients( tEvaluationPoint );

                            Matrix< DDRMat > tForwardObjectiveGradient  = this->evaluate_objective_gradients( tEvaluationName + " fd_forward" );
                            Matrix< DDRMat > tForwardConstraintGradient = this->evaluate_constraint_gradients( tEvaluationName + " fd_forward" );

                            // Backward
                            this->set_sensitivity_analysis_type( SA_Type::backward );
                            this->compute_design_criteria_gradients( tEvaluationPoint );

                            Matrix< DDRMat > tBackwardObjectiveGradient  = this->evaluate_objective_gradients( tEvaluationName + " fd_backward" );
                            Matrix< DDRMat > tBackwardConstraintGradient = this->evaluate_constraint_gradients( tEvaluationName + " fd_backward" );

                            // Central
                            this->output_variables( ( tForwardObjectiveGradient + tBackwardObjectiveGradient ) / 2,
                                    "objective_gradients" + tEvaluationName + " fd_central" );
                            this->output_variables( ( tForwardConstraintGradient + tBackwardConstraintGradient ) / 2,
                                    "constraint_gradients" + tEvaluationName + " fd_central" );
                        }
                    }
                }
            }

            // Close file
            close_hdf5_file( mFileID );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_Sweep::output_objectives_constraints( std::string aEvaluationName )
        {
            // Output
            if ( mEvaluateObjectives )
            {
                const Matrix< DDRMat >& tObjectives = this->get_objectives();
                this->output_variables( tObjectives, "objectives" + aEvaluationName );
            }

            if ( mEvaluateConstraints )
            {
                const Matrix< DDRMat >& tConstraints = this->get_constraints();
                this->output_variables( tConstraints, "constraints" + aEvaluationName );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Algorithm_Sweep::evaluate_objective_gradients( std::string aEvaluationName )
        {
            if ( mEvaluateObjectiveGradients )
            {
                // Calculate the objective gradients; note: a copy is created as the gradient
                // vector is altered below
                Matrix< DDRMat > tObjectiveGradients = this->get_objective_gradients();

                // Extract sensitivities with respect to advs specified for FD
                if ( tObjectiveGradients.numel() > mFiniteDifferenceADVs.numel() )
                {
                    // Allocate temporary matrix (needed as adv indices in mFiniteDifferenceADVs may not be sorted
                    Matrix< DDRMat > tObjectiveGradientsExtracted( mFiniteDifferenceADVs.numel(), 1 );

                    // Extract sensitivities of requested ADVs
                    for ( uint tIndex = 0; tIndex < mFiniteDifferenceADVs.numel(); ++tIndex )
                    {
                        // get ADV index
                        uint tAdvIndex = mFiniteDifferenceADVs( tIndex );

                        MORIS_ASSERT( tAdvIndex < tObjectiveGradients.numel(),
                                "Algorithm_Sweep::evaluate_objective_gradients - requested adv index does not exist.\n" );

                        tObjectiveGradientsExtracted( tIndex ) = tObjectiveGradients( tAdvIndex );
                    }

                    // replace full derivative matrix with temporary one
                    tObjectiveGradients = tObjectiveGradientsExtracted;
                }

                // Output
                this->output_variables( tObjectiveGradients, "objective_gradients" + aEvaluationName );

                // Return
                return tObjectiveGradients;
            }
            return { {} };
        }

        //--------------------------------------------------------------------------------------------------------------

        Matrix< DDRMat >
        Algorithm_Sweep::evaluate_constraint_gradients( std::string aEvaluationName )
        {
            if ( mEvaluateConstraintGradients )
            {
                // Calculate the constraint gradients; note: a copy is created as the gradient
                // matrix is altered below
                Matrix< DDRMat > tConstraintGradients = this->get_constraint_gradients();

                // Extract sensitivities with respect to advs specified for FD
                if ( tConstraintGradients.n_cols() > mFiniteDifferenceADVs.numel() )
                {
                    // Number of constraints
                    uint tNumberOfConstraints = tConstraintGradients.n_rows();

                    // Allocate temporary matrix (needed as adv indices in mFiniteDifferenceADVs may not be sorted
                    Matrix< DDRMat > tConstraintGradientsExtracted( tNumberOfConstraints, mFiniteDifferenceADVs.numel(), 1 );

                    // Extract sensitivities of requested ADVs
                    for ( uint tIndex = 0; tIndex < mFiniteDifferenceADVs.numel(); ++tIndex )
                    {
                        for ( uint tConst = 0; tConst < tNumberOfConstraints; ++tConst )
                        {
                            // get ADV index
                            uint tAdvIndex = mFiniteDifferenceADVs( tIndex );

                            MORIS_ASSERT( tAdvIndex < tConstraintGradients.n_cols(),
                                    "Algorithm_Sweep::evaluate_constraint_gradients - requested adv index does not exist.\n" );

                            tConstraintGradientsExtracted( tConst, tIndex ) = tConstraintGradients( tConst, tAdvIndex );
                        }
                    }

                    // replace full derivative matrix with temporary one
                    tConstraintGradients = tConstraintGradientsExtracted;
                }

                // Output
                this->output_variables( tConstraintGradients, "constraint_gradients" + aEvaluationName );

                // Return
                return tConstraintGradients;
            }
            return { {} };
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_Sweep::output_variables(
                Matrix< DDRMat > aVariables,
                std::string      aFullEvaluationName )
        {
            // Write status
            herr_t tStatus = 0;

            // Save
            if ( mSave )
            {
                moris::save_matrix_to_hdf5_file( mFileID, aFullEvaluationName, aVariables, tStatus );
            }

            // Print
            if ( mPrint )
            {
                moris::print( aVariables, aFullEvaluationName );
            }
        }

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace opt
}    // namespace moris
