/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_LBFGS.cpp
 *
 */

#include "cl_OPT_Algorithm_LBFGS.hpp"
#include "cl_Communication_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "fn_Parsing_Tools.hpp"
#include "fn_norm.hpp"
#include "cl_Fortran.hpp"

#ifdef FORT_NO_
#define _FORTRAN( a ) a
#elif F77ADD_
#define _FORTRAN( a ) a##_
#elif F77ADD__
#define _FORTRAN( a ) a##__
#endif

#ifdef MORIS_HAVE_LBFGS
extern "C" {
// L-BFGS-B (a Fortran implementation of BFGS) function declaration
void setulb_(
        int&    n,
        int&    m,
        double* x,
        double* l,
        double* u,
        int*    nbd,
        double& f,
        double* g,
        double& mNormDrop,
        double& mGradTolerance,
        double* wa,
        int*    iwa,
        char*   task,
        int&    iprint,
        char*   csave,
        int*    lsave,
        int*    isave,
        double* dsave );
}
#endif

//----------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_LBFGS::Algorithm_LBFGS( ParameterList aParameterList )
                : mMaxIt( aParameterList.get< sint >( "max_its" ) )
                , mLBFGSprint( aParameterList.get< sint >( "internal_lbfgs_print_severity" ) )
                , mNumCorrections( aParameterList.get< sint >( "num_corr" ) )
                , mNumberOfFunctionEvals( aParameterList.get< sint >( "num_function_evaluations" ) )
                , mNormDrop( aParameterList.get< real >( "norm_drop" ) )
                , mGradTolerance( aParameterList.get< real >( "grad_tol" ) )
        {
#ifndef MORIS_HAVE_LBFGS
            MORIS_ERROR( false, "MORIS was compiled without LBFGS support" );
#endif

            // convert input parameters ti matrix
            string_to_mat< DDRMat >( aParameterList.get< std::string >( "step_size" ), mStepSize );
            string_to_mat< DDUMat >( aParameterList.get< std::string >( "outer_iteration_index" ), mOuterIterationIndex );
            string_to_mat< DDUMat >( aParameterList.get< std::string >( "number_inner_iterations" ), mNumberOfInnerIterations );

            // check if the input matrices have the same size
            MORIS_ASSERT( mStepSize.numel() == mOuterIterationIndex.numel() and mOuterIterationIndex.numel() == mNumberOfInnerIterations.numel(),
                    "Algorithm_LBFGS::Algorithm_LBFGS failed, three inputs step_size,step_size_index and number_inner_iterations must have the same size" );
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_LBFGS::~Algorithm_LBFGS()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        uint
        Algorithm_LBFGS::solve(
                uint                       aCurrentOptAlgInd,
                std::shared_ptr< Problem > aOptProb )
        {
            // Trace optimization
            Tracer tTracer( "OptimizationAlgorithm", "LBFGS", "Solve" );

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
                // Run gcmma algorithm
                this->lbfgs_solve();

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
        Algorithm_LBFGS::lbfgs_solve()
        {
#ifdef MORIS_HAVE_LBFGS
            // get number of design variables
            int tNumAdvs = mProblem->get_num_advs();

            // NOTE:  this is done in order to be able to use native functions on the input data that will be passed on to Fortran
            // algorithm inputs in terms of cell
            moris::Vector< int >    nbdCell( tNumAdvs, 2 );    // This algorithm parameter is hard-coded to assume that the variables are bounded
            moris::Vector< double > gCell( tNumAdvs );
            moris::Vector< double > waCell( ( 2 * mNumCorrections + 5 ) * tNumAdvs + 11 * mNumCorrections * mNumCorrections + 8 * mNumCorrections );
            moris::Vector< int >    iwaCell( 3 * tNumAdvs );
            moris::Vector< int >    isaveCell( 44 );
            moris::Vector< double > dsaveCell( 29 );

            // initialize and allocate memory for algorithm inputs/outputs
            int*    nbd = nbdCell.memptr();
            double  f   = 0;
            double* g   = gCell.memptr();
            double* wa  = waCell.memptr();
            int*    iwa = iwaCell.memptr();

            int iprint = mLBFGSprint;    // Prevents the algorithm from printing anything

            moris::fortran::LOGICAL lsave[ 4 ] = { false };

            int*    isave = isaveCell.memptr();
            double* dsave = dsaveCell.memptr();

            // convert user-defind required norm drop to LBFGS specific definition
            real tLbfgsNormDrop = mNormDrop / MORIS_REAL_EPS;

            // starts the algorithm, if initialized static it must  be of size 61
            char task[ 60 + 1 ];
            task[ 60 ] = '\0';

            // initialize empty str
            char csave[ 60 + 1 ];
            csave[ 60 ] = '\0';

            // get optimization variables
            Matrix< DDRMat >& tAdvs = mProblem->get_advs();

            // upper and lower bounds (fixed and prescribed)
            const Matrix< DDRMat >& tLowerBounds = mProblem->get_lower_bounds();
            const Matrix< DDRMat >& tUpperBounds = mProblem->get_upper_bounds();

            // check for consistent length of variable and bound vectors
            MORIS_ERROR( tAdvs.numel() == tLowerBounds.numel() && tAdvs.numel() == tUpperBounds.numel(),
                    "Algorithm_LBFGS::lbfgs_solve - inconsistent length of ADV and lower/upper bound vectors." );

            // create arrays for temporary upper and lower bounds
            Matrix< DDRMat > tTmpLowerBounds( tNumAdvs, 1 );
            Matrix< DDRMat > tTmpUpperBounds( tNumAdvs, 1 );

            // step size adjustment and counter to determine the which step size to use
            uint        iCounter                 = 0;
            moris::real tStep                    = mStepSize( iCounter );
            moris::uint tNumberOfInnerIterations = mNumberOfInnerIterations( iCounter );

            // indicator that problem has converged
            bool tIsConverged = false;

            // outer iteration for lbfgs
            for ( uint iOuterIteration = 0; iOuterIteration < (uint)mMaxIt; iOuterIteration++ )
            {
                // adjust the step size and number of inner iterations
                if ( iCounter + 1 < mStepSize.numel() )
                {
                    if ( iOuterIteration == mOuterIterationIndex( iCounter + 1 ) )
                    {
                        // adjust step size
                        tStep = mStepSize( iCounter + 1 );

                        // adjust number of inner iterations
                        tNumberOfInnerIterations = mNumberOfInnerIterations( iCounter + 1 );

                        // increase counter
                        iCounter++;
                    }
                }

                MORIS_LOG_INFO( "Outer Iteration: %d  StepSize = %e  Num Inner Iterations = %d",
                        iOuterIteration + 1,
                        tStep,
                        tNumberOfInnerIterations );

                // determine user defined upper and lower bounds
                for ( uint iADV = 0; iADV < (uint)tNumAdvs; iADV++ )
                {
                    tAdvs( iADV ) = std::max( tLowerBounds( iADV ), std::min( tAdvs( iADV ), tUpperBounds( iADV ) ) );

                    tTmpLowerBounds( iADV ) = std::max( tLowerBounds( iADV ), tAdvs( iADV ) - tStep * ( tUpperBounds( iADV ) - tLowerBounds( iADV ) ) );
                    tTmpUpperBounds( iADV ) = std::min( tUpperBounds( iADV ), tAdvs( iADV ) + tStep * ( tUpperBounds( iADV ) - tLowerBounds( iADV ) ) );
                }

                // extract pointers for optimization variables and temporary upper and lower bounds
                double* x = tAdvs.data();
                double* l = tTmpLowerBounds.data();
                double* u = tTmpUpperBounds.data();

                // restart the optimization in every outer iter
                strcpy( task, "START" );
                strcpy( csave, "" );

                // convert to fortran string (size is one less than allocated)
                moris::fortran::CHARACTER tTaskFortran( task, 60 );
                moris::fortran::CHARACTER tCsaveFortran( csave, 60 );

                // log total and inner iterations
                MORIS_LOG_INFO( "Total iteration: %d   Inner iteration: %d", mOptIter + 1, 1 );

                // set convergence flag to false before inner iterations
                tIsConverged = false;

                // inner loop to go through internal iterations
                uint iInnerIteration = 0;
                while ( iInnerIteration < tNumberOfInnerIterations )
                {
                    // call the Fortran subroutine
                    setulb_(
                            tNumAdvs,
                            mNumCorrections,
                            x,
                            l,
                            u,
                            nbd,
                            f,
                            g,
                            tLbfgsNormDrop,
                            mGradTolerance,
                            wa,
                            iwa,
                            tTaskFortran,
                            iprint,
                            tCsaveFortran,
                            lsave,
                            isave,
                            dsave );

                    // evaluate task string
                    if ( strncmp( task, "ERROR", 5 ) == 0 )
                    {
                        MORIS_ERROR( false,
                                "Algorithm_LBFGS::lbfgs_solve - %s",
                                task );
                    }
                    if ( strncmp( task, "CONV", 4 ) == 0 )
                    {
                        // set convergence flag to true
                        tIsConverged = true;

                        // log task
                        MORIS_LOG_INFO( "Task: %s", task );

                        // leave inner iteration
                        break;
                    }
                    else if ( strncmp( task, "NEW_X", 5 ) == 0 )
                    {
                        // increase inner iteration count
                        iInnerIteration++;

                        // increase total iteration count
                        mOptIter++;

                        // log total and inner iterations
                        MORIS_LOG_INFO( "Total iterations: %d   Inner iteration: %d  ",
                                mOptIter + 1,
                                iInnerIteration + 1 );

                        // log number of function calls, objective and projected gradient norm
                        MORIS_LOG_INFO( "Function Evals = %d   Objective = %e   |proj g| = %e",
                                isave[ 33 ],
                                f,
                                dsave[ 12 ] );

                        // check whether stop criteria met

                        // terminate if the total number of f and g evaluations exceeds limit
                        if ( isave[ 33 ] >= mNumberOfFunctionEvals )
                        {
                            MORIS_LOG_INFO( "Total number of function calls has been reached" );

                            break;
                        }

                        // terminate if  |proj g|/(1+|f|) < 1.0d-10, where "proj g" denoted the projected gradient
                        if ( dsave[ 12 ] <= 1e-10 * ( 1.0 + std::abs( f ) ) )
                        {
                            // set convergence flag to true
                            tIsConverged = true;

                            MORIS_LOG_INFO( "The projected gradient is sufficiently small" );

                            break;
                        }
                    }
                    else if ( strncmp( task, "FG", 2 ) == 0 )
                    {
                        // compute objective
                        this->func( mOptIter, x, f );

                        // compute gradients
                        this->grad( x, g );
                    }
                }

                // check if the LBFGS has converged
                if ( tIsConverged )
                {
                    // check if ADVs are within temporary bounds or on prescribed bounds
                    for ( uint iADV = 0; iADV < (uint)tNumAdvs; iADV++ )
                    {
                        // set convergence flag for current ADV to false
                        bool tAdvConverged = false;

                        // option 1: design variable is within the user defined bounds
                        if (
                                ( x[ iADV ] - tTmpLowerBounds( iADV ) ) >= MORIS_REAL_EPS and    //
                                ( tTmpUpperBounds( iADV ) - x[ iADV ] ) >= MORIS_REAL_EPS )
                        {
                            tAdvConverged = true;
                        }

                        // option 2: design variable is on prescribed upper bound
                        else if ( std::abs( x[ iADV ] - tUpperBounds( iADV ) ) < MORIS_REAL_EPS )
                        {
                            tAdvConverged = true;
                        }

                        // option 3: design variable is on prescribed lower bound
                        else if ( std::abs( x[ iADV ] - tLowerBounds( iADV ) ) < MORIS_REAL_EPS )
                        {
                            tAdvConverged = true;
                        }

                        // check if current ADV is converged
                        if ( !tAdvConverged )
                        {
                            MORIS_LOG_INFO( "Solution of inner iterations not within or on upper/lower bounds" );

                            tIsConverged = false;
                            break;
                        }
                    }

                    // if all the adv satisfy the bounds then we have converged and outer loop breaks
                    if ( tIsConverged == true )
                    {
                        MORIS_LOG_INFO( "LBFGS has converged" );
                        break;
                    }
                }
            }
#endif
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_LBFGS::func( int aIter, double* aAdv, double& aObjval )
        {
            // Update the ADV matrix
            Matrix< DDRMat > tADVs( aAdv, mProblem->get_num_advs(), 1, false, true );

            // Write restart file
            this->write_advs_to_file( tADVs );

            // Recruit help from other procs and solve for criteria
            this->compute_design_criteria( tADVs );

            // Convert outputs from type MORIS
            aObjval = mProblem->get_objectives()( 0 );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_LBFGS::grad(
                double* aAdv,
                double* aD_Obj )
        {
            // Update the ADV matrix
            Matrix< DDRMat > tADVs( aAdv, mProblem->get_num_advs(), 1, false, true );

            // Compute design criteria gradients
            this->compute_design_criteria_gradients( tADVs );

            // copy the gradients into the memebr data in problem
            auto tD_Obj = mProblem->get_objective_gradients().data();
            std::copy( tD_Obj, tD_Obj + mProblem->get_num_advs(), aD_Obj );
        }

        //--------------------------------------------------------------------------------------------------------------

        void
        Algorithm_LBFGS::printresult()
        {
            std::fprintf( stdout, " \nResult of LBFGS\n" );
        }

        //--------------------------------------------------------------------------------------------------------------
    }    // namespace opt
}    // namespace moris
