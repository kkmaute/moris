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
        int*    n,
        int*    m,
        double* x,
        double* l,
        double* u,
        int*    nbd,
        double* f,
        double* g,
        double* mNormDrop,
        double* mGradTolerance,
        double* wa,
        int*    iwa,
        char*   task,
        int*    iprint,
        char*   csave,
        bool*   lsave,
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
            int n = mProblem->get_num_advs();

            // Note that these pointers are deleted by the the Arma and Eigen
            // libraries themselves.
            auto x = mProblem->get_advs().data();
            auto l = mProblem->get_lower_bounds().data();
            auto u = mProblem->get_upper_bounds().data();

            // NOTE:  this is done in order to be able to use native functions on the input data that will be passed on to Fortran
            // algorithm inputs in terms of cell
            moris::Cell< int >    nbdCell( n, 2 );    // This algorithm parameter is hard-coded to assume that the variables are bounded
            moris::Cell< double > gCell( n );
            moris::Cell< double > waCell( ( 2 * mNumCorrections + 5 ) * n + 11 * mNumCorrections * mNumCorrections + 8 * mNumCorrections );
            moris::Cell< int >    iwaCell( 3 * n );
            moris::Cell< int >    isaveCell( 44 );
            moris::Cell< double > dsaveCell( 29 );

            // initialize and allocate memory for algorithm inputs/outputs
            int*    nbd = nbdCell.memptr();
            double  f   = 0;
            double* g   = gCell.memptr();
            double* wa  = waCell.memptr();
            int*    iwa = iwaCell.memptr();

            int iprint = mLBFGSprint;    // Prevents the algorithm from printing anything

            char csave[ 61 ];
            bool lsave[ 4 ] = { false };

            int*    isave = isaveCell.memptr();
            double* dsave = dsaveCell.memptr();

            // starts the algorithm , if initialized static it must  be of size 61
            char* task = new char[ 61 ]{ 'S', 'T', 'A', 'R', 'T', '\0' };

            // pad the string with the spaces as it is the convention in fortran
            std::fill( &task[ 0 ] + std::strlen( &task[ 0 ] ), &task[ 0 ] + 60, ' ' );

            // upper and lower bounds user defined will be overwritten later
            Matrix< DDRMat > tLowerBoundsUserDefined = mProblem->get_lower_bounds();
            Matrix< DDRMat > tUpperBoundsUserDefined = mProblem->get_upper_bounds();

            // upper and lower bounds fixed and prescribed
            Matrix< DDRMat > tLowerBounds = mProblem->get_lower_bounds();
            Matrix< DDRMat > tUpperBounds = mProblem->get_upper_bounds();

            // step size adjustment and counter to determine the which step size to use
            uint        iCounter                 = 0;
            moris::real tStep                    = mStepSize( iCounter );
            moris::uint tNumberOfInnerIterations = mNumberOfInnerIterations( iCounter );

            // indicator that problem has converged
            bool tIsConverged = false;

            // outer iteration for lbfgs
            for ( size_t iOuterIteration = 0; iOuterIteration < (uint)mMaxIt; iOuterIteration++ )
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

                // determine user defined upper and lower bounds
                for ( size_t iADV = 0; iADV < tUpperBoundsUserDefined.numel(); iADV++ )
                {
                    tLowerBoundsUserDefined( iADV ) = std::max( tLowerBounds( iADV ), x[ iADV ] - tStep * ( tUpperBounds( iADV ) - tLowerBounds( iADV ) ) );
                    tUpperBoundsUserDefined( iADV ) = std::min( tUpperBounds( iADV ), x[ iADV ] + tStep * ( tUpperBounds( iADV ) - tLowerBounds( iADV ) ) );
                }

                // extract pointers
                l = tLowerBoundsUserDefined.data();
                u = tUpperBoundsUserDefined.data();

                // value to store previous optimization iteration
                double f_previous = 0.0;

                // inner loop to go through internal iterations
                for ( size_t iInnerIteration = 0; iInnerIteration < tNumberOfInnerIterations; iInnerIteration++ )
                {
                    // call the Fortran subroutine
                    setulb_( &n,
                            &mNumCorrections,
                            x,
                            l,
                            u,
                            nbd,
                            &f,
                            g,
                            &mNormDrop,
                            &mGradTolerance,
                            wa,
                            iwa,
                            task,
                            &iprint,
                            csave,
                            lsave,
                            isave,
                            dsave );

                    // only assign f_previous at iteration 1
                    f_previous = iInnerIteration == 1 ? f : f_previous;

                    if ( strncmp( task, "FG", 2 ) == 0 || strncmp( task, "NEW_X", 5 ) == 0 )
                    {
                        // call to compute objective
                        this->func( mOptIter, x, f );

                        // call to compute gradients
                        this->grad( x, g );
                    }

                    // increment the iteration
                    mOptIter++;
                }

                // check if the LBFGS is converging through objective values
                if ( std::abs( f - f_previous ) / std::max( { f, f_previous, 1.0 } ) < mNormDrop * MORIS_REAL_EPS )
                {
                    // check if the bounds are satisfied
                    for ( size_t iADV = 0; iADV < tUpperBounds.numel(); iADV++ )
                    {
                        // option 1: the design variables are within the user defined bounds
                        if ( x[ iADV ] > tLowerBoundsUserDefined( iADV ) and x[ iADV ] < tUpperBoundsUserDefined( iADV ) )
                        {
                            tIsConverged = true;
                            continue;
                        }

                        // the design variable hits the the user defined upper bound then the use defined bound and problem bound should be the same
                        else if ( std::abs( x[ iADV ] - tUpperBoundsUserDefined( iADV ) ) < MORIS_REAL_EPS )
                        {
                            if ( std::abs( tUpperBounds( iADV ) - tUpperBoundsUserDefined( iADV ) ) < MORIS_REAL_EPS )
                            {
                                tIsConverged = true;
                                continue;
                            }

                            // break if this is the case
                            tIsConverged = false;
                            break;
                        }

                        // the design variable hits the the user defined lower bound then the use defined bound and problem bound should be the same
                        else if ( std::abs( x[ iADV ] - tLowerBoundsUserDefined( iADV ) ) < MORIS_REAL_EPS )
                        {
                            if ( std::abs( tLowerBounds( iADV ) - tLowerBoundsUserDefined( iADV ) ) < MORIS_REAL_EPS )
                            {
                                tIsConverged = true;
                                continue;
                            }

                            // break if this is the case
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

            // delete the pointer
            delete[] task;
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
