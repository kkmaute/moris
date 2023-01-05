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

//----------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_LBFGS::Algorithm_LBFGS( ParameterList aParameterList )
                : mMaxIt( aParameterList.get< sint >( "max_its" ) )
                , mNumCorrections( aParameterList.get< sint >( "num_corr" ) )
                , mNormDrop( aParameterList.get< real >( "norm_drop" ) )
                , mGradTolerance( aParameterList.get< real >( "grad_tol" ) )
                , mLBFGSprint( aParameterList.get< sint >( "internal_lbfgs_print_severity" ) )
        {
            // convert input parameters ti matrix
            string_to_mat< DDRMat >( aParameterList.get< std::string >( "step_size" ), mStepSize );
            string_to_mat< DDUMat >( aParameterList.get< std::string >( "step_size_index" ), mStepSizeIndex );
            string_to_mat< DDUMat >( aParameterList.get< std::string >( "number_inner_iterations" ), mNumberOfInnerIterations );

            // check if the input matrices have the same size
            MORIS_ASSERT( mStepSize.numel() == mStepSizeIndex.numel() and mStepSizeIndex.numel() == mNumberOfInnerIterations.numel(),
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
            // get number of design variables
            int tNumAdvs = mProblem->get_num_advs();

            // NOTE:  this is done in order to be able to use native functions on the input data that will be passed on to Fortran
            // algorithm inputs in terms of cell
            moris::Cell< int >    nbdCell( tNumAdvs, 2 );    // This algorithm parameter is hard-coded to assume that the variables are bounded
            moris::Cell< double > gCell( tNumAdvs );
            moris::Cell< double > waCell( ( 2 * mNumCorrections + 5 ) * tNumAdvs + 11 * mNumCorrections * mNumCorrections + 8 * mNumCorrections );
            moris::Cell< int >    iwaCell( 3 * tNumAdvs );
            moris::Cell< int >    isaveCell( 44 );
            moris::Cell< double > dsaveCell( 29 );

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
                // adjust the step size if there is more than 1 step size through optimization
                // NOTE: step size only should be adjusted in outer iterations
                if ( iCounter + 1 < mStepSize.numel() )
                {
                    if ( mOptIter == mStepSizeIndex( iCounter + 1 ) )
                    {
                        // adjust the step size
                        tStep                    = mStepSize( iCounter + 1 );
                        tNumberOfInnerIterations = mNumberOfInnerIterations( iCounter + 1 );
                        iCounter++;
                    }
                }

                MORIS_LOG_INFO( "%d Outer Iteration  StepSize = %e  Num Inner Iterations = %d",
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

                // value to store previous optimization iteration
                double f_previous = 0.0;

                // restart the optimization in every outer iter
                strcpy( task, "START" );
                strcpy( csave, "" );

                // convert to fortran string (size is one less than allocated)
                moris::fortran::CHARACTER tTaskFortran( task, 60 );
                moris::fortran::CHARACTER tCsaveFortran( csave, 60 );

                // inner loop to go through internal iterations
                // for ( uint iInnerIteration = 0; iInnerIteration < tNumberOfInnerIterations; iInnerIteration++ )
                uint iInnerIteration = 0;
                while ( iInnerIteration < tNumberOfInnerIterations )
                {
                    std::cout<<"xxx task: "<<std::to_string(iInnerIteration)<< "  : " <<task<<std::endl;
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
                            mNormDrop,
                            mGradTolerance,
                            wa,
                            iwa,
                            tTaskFortran,
                            iprint,
                            tCsaveFortran,
                            lsave,
                            isave,
                            dsave );

                    // only assign f_previous at iteration 1
                    f_previous = iInnerIteration == 0 ? f : f_previous;

                    // evaluate task string
                    if ( strncmp( task, "ERROR", 5 ) == 0 )
                    {
                        MORIS_ERROR( false,
                                "Algorithm_LBFGS::lbfgs_solve - %s",
                                task );
                    }
                    if ( strncmp( task, "CONV", 4 ) == 0 )
                    {
                        MORIS_LOG_INFO( "%s", task );

                        // leave inner iteration
                        break;
                    }
                    else if ( strncmp( task, "NEW_X", 5 ) == 0 )
                    {
                        iInnerIteration++;
                        // do nothing for right now; could check for convergence (see driver2.f)
                    }
                    else if ( strncmp( task, "FG", 2 ) == 0 )
                    {
                        // call to compute objective
                        this->func( mOptIter, x, f );

                        // call to compute gradients
                        this->grad( x, g );
                    }

                    // increment the iteration
                    mOptIter++;
                }

                // compute reference value for convergence check
                real tReferenceObjective = std::max( { std::abs( f ), std::abs( f_previous ), 1.0 } );

                // check if the LBFGS has converged wrt objective values
                if ( std::abs( f - f_previous ) / tReferenceObjective < mNormDrop * MORIS_REAL_EPS )
                {
                    // set convergence to true
                    tIsConverged = true;

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
