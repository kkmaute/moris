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

            gLogger.set_iteration( "OptimizationManager", LOGGER_NON_SPECIFIC_ENTITY_TYPE, "Perform", tOptIter );

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

            // run the algorithm until converges
            while ( ( strncmp( task, "FG", 2 ) == 0 ) ||       //
                    ( strncmp( task, "NEW_X", 5 ) == 0 ) ||    //
                    ( strncmp( task, "START", 5 ) == 0 ) )
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

                if ( strncmp( task, "FG", 2 ) == 0 )
                {
                    // call to compute objective
                    this->func( mOptIter, x, f );

                    // call to compute gradients
                    this->grad( x, g );
                }

                // one iteration of algorithm has concluded
                if ( strncmp( task, "NEW_X", 5 ) == 0 )
                {
                    // set optimization iteration counter
                    mOptIter = isave[ 29 ];
                }

                // exit loop if maximum iterations have been achieved
                if ( isave[ 29 ] == mMaxIt )
                {
                    break;
                }
            }

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
