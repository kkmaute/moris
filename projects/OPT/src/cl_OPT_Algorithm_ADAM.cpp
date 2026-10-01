/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_ADAM.cpp
 *
 */

#include "cl_OPT_Algorithm_ADAM.hpp"
#include "cl_Communication_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_norm.hpp"
#include "cl_Fortran.hpp"

//----------------------------------------------------------------------------------------------------------------------

namespace moris::opt
{

    //--------------------------------------------------------------------------------------------------------------

    Algorithm_ADAM::Algorithm_ADAM( const Parameter_List& aParameterList )
            : Algorithm()
            , mMaxIt( aParameterList.get< sint >( "max_its" ) )
            , mNormDrop( aParameterList.get< real >( "norm_drop" ) )
            , mGradTolerance( aParameterList.get< real >( "grad_tol" ) )
            , mMomentumBeta1( aParameterList.get< real >( "gradient_decay_factor" ) )
            , mMomentumBeta2( aParameterList.get< real >( "squared_gradient_decay_factor" ) )
            , mLearningRate( aParameterList.get< real >( "step_size" ) )
            , mEpsilon( aParameterList.get< real >( "epsilon" ) )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Algorithm_ADAM::~Algorithm_ADAM()
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint
    Algorithm_ADAM::solve(
            uint                       aCurrentOptAlgInd,
            std::shared_ptr< Problem > aOptProb )
    {
        // Trace optimization
        Tracer tTracer( "OptimizationAlgorithm", "ADAM", "Solve" );

        // running status has to be wait when starting a solve
        mRunning = opt::Task::wait;

        mCurrentOptAlgInd = aCurrentOptAlgInd;    // set index of current optimization algorithm
        mProblem          = aOptProb;             // set the member variable mProblem to aOptProb

        // Set optimization iteration index for restart
        if ( mRestartIndex > 0 )
        {
            gLogger.set_opt_iteration( mRestartIndex );
        }

        // Solve optimization problem
        if ( par_rank() == 0 )
        {
            // Run ADAM algorithm
            this->ADAM_solve();

            // Communicate that optimization has finished
            mRunning = opt::Task::exit;

            this->communicate_running_status();
        }
        else
        {
            // Run dummy solve
            this->dummy_solve();
        }

        uint tOptIter = gLogger.get_opt_iteration();

        gLogger.set_iteration( "OPT", "Manager", "Perform", tOptIter );

        return tOptIter;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ADAM::ADAM_solve()
    {
        // get optimization variables and their bounds
        Vector< real >& tADVs    = mProblem->get_advs();
        uint            tNumADVs = tADVs.size();

        Vector< real >& tLowerBounds = mProblem->get_lower_bounds();
        Vector< real >& tUpperBounds = mProblem->get_upper_bounds();

        // Initialize previous gradient and squared gradient for momentum calculation
        Vector< real > tPrevFirstMoment( tNumADVs, 0.0 );     // Previous gradient for momentum calculation
        Vector< real > tPrevSecondMoment( tNumADVs, 0.0 );    // Previous squared gradient for momentum calculation

        // Initialize first and second moments
        real           tFirstMoment;
        real           tSecondMoment;
        real           tBiasCorrectedFirstMoment;
        real           tBiasCorrectedSecondMoment;
        real           tGradient;
        real           tBeta1Power = 1.0;
        real           tBeta2Power = 1.0;
        Vector< real > tDeltaADVs( tNumADVs );    // Change in ADVs for convergence check

        // Loop over maximum number of iterations
        for ( uint iIt = 0; iIt < mMaxIt; ++iIt )
        {
            // Recruit help from other procs and solve for criteria
            this->compute_design_criteria( tADVs );

            // Compute design criteria gradients
            this->compute_design_criteria_gradients( tADVs );

            // Compute loss gradient
            const Matrix< DDRMat >& tGradients = this->get_objective_gradients();

            // Check for convergence based on norm of the gradient
            if ( norm( tGradients ) < mGradTolerance )
            {
                MORIS_LOG_INFO( "Convergence achieved based on gradient tolerance at iteration %u.", iIt + 1 );
                return;
            }

            // Update the beta powers
            tBeta1Power *= mMomentumBeta1;
            tBeta2Power *= mMomentumBeta2;

            // Update ADVs based on ADAM update rule
            for ( uint iDV = 0; iDV < tADVs.size(); ++iDV )
            {
                // Get the gradient
                tGradient = tGradients( iDV );

                // Compute first and second moments

                tFirstMoment  = mMomentumBeta1 * tPrevFirstMoment( iDV ) + ( 1 - mMomentumBeta1 ) * tGradient;
                tSecondMoment = mMomentumBeta2 * tPrevSecondMoment( iDV ) + ( 1 - mMomentumBeta2 ) * tGradient * tGradient;

                // Bias correction for both moments
                tBiasCorrectedFirstMoment  = tFirstMoment / ( 1.0 - tBeta1Power );
                tBiasCorrectedSecondMoment = tSecondMoment / ( 1.0 - tBeta2Power );

                // Store the delta in ADVs
                tDeltaADVs( iDV ) = mLearningRate * tBiasCorrectedFirstMoment / ( std::sqrt( tBiasCorrectedSecondMoment ) + mEpsilon );

                // Store current gradient and squared gradient for next iteration
                tPrevFirstMoment( iDV )  = tFirstMoment;
                tPrevSecondMoment( iDV ) = tSecondMoment;
            }

            // Update ADVs with box constraints
            for ( uint iDV = 0; iDV < tNumADVs; ++iDV )
            {
                real tOldADV = tADVs( iDV );

                // Proposed ADAM update
                real tProposedADV = tOldADV - tDeltaADVs( iDV );

                // Project onto box constraints
                tADVs( iDV ) = std::max( tLowerBounds( iDV ), std::min( tUpperBounds( iDV ), tProposedADV ) );

                // Store actual ADV change
                tDeltaADVs( iDV ) = tOldADV - tADVs( iDV );
            }

            // Write ADVs to file
            this->write_advs_to_file( tADVs );

            PRINT( mProblem->get_objectives() );

            // Check for convergence based on the change in ADVs
            if ( std::sqrt( std::accumulate( tDeltaADVs.begin(), tDeltaADVs.end(), 0.0, []( real a, real b ) { return a + b * b; } ) ) < mNormDrop )
            {
                MORIS_LOG_INFO( "Convergence achieved based on step size tolerance at iteration %u.", iIt + 1 );
                return;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ADAM::printresult()
    {
        std::fprintf( stdout, " \nResult of ADAM\n" );
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::opt
