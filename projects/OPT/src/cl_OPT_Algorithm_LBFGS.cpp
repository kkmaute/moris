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

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#ifdef FORT_NO_
#define _FORTRAN(a) a
#elif  F77ADD_
#define _FORTRAN(a) a ## _
#elif F77ADD__
#define _FORTRAN(a) a ## __
#endif

extern "C"
{
    // L-BFGS-B (a Fortran implementation of BFGS) function declaration
    void setulb_(
            int* n, int* m, double* x, double* l, double* u, int* nbd,
            double* f, double* g, double* mNormDrop, double* mGradTolerance,
            double* wa, int* iwa, char* task, int* iprint, char* csave,
            bool* lsave, int* isave, double* dsave);
}

//----------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_LBFGS::Algorithm_LBFGS(ParameterList aParameterList)
                : mMaxIt(aParameterList.get< sint >( "max_its" )),
                  mNumCorrections(aParameterList.get< sint >( "num_corr" )),
                  mNormDrop(aParameterList.get< real >( "norm_drop" )),
                  mGradTolerance(aParameterList.get< real >( "grad_tol" ))
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_LBFGS::~Algorithm_LBFGS()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Algorithm_LBFGS::solve( uint aCurrentOptAlgInd, std::shared_ptr<Problem> aOptProb )
        {
            // Trace optimization
            Tracer tTracer( "OptimizationAlgorithm", "LBFGS", "Solve" );

            mCurrentOptAlgInd = aCurrentOptAlgInd;  // set index of current optimization algorithm
            mProblem          = aOptProb;           // set the member variable mProblem to aOptProb

            // Set optimization iteration index for restart
            if ( mRestartIndex > 0)
            {
                gLogger.set_opt_iteration( mRestartIndex );
            }

            int n = mProblem->get_num_advs(); // number of design variables

            // Note that these pointers are deleted by the the Arma and Eigen
            // libraries themselves.
            auto x = mProblem->get_advs().data();
            auto l = mProblem->get_lower_bounds().data();
            auto u = mProblem->get_upper_bounds().data();

            // This algorithm parameter is hard-coded to assume that the variables are bounded
            int* nbd = new int[n];
            std::fill( nbd, nbd + n, 2 );

            // initialize and allocate memory for algorithm inputs/outputs
            double f       = 0;
            double* g      = new double[n];
            double* wa     = new double[(2*mNumCorrections + 5)*n + 11*mNumCorrections*mNumCorrections + 8*mNumCorrections];
            int* iwa       = new int[3 * n];
            int iprint     = -1;      // Prevents the algorithm from printing anything
            char task[61]  = "START"; // starts the algorithm
            char fg_start[61] = "FG_START                                                    "; // needed for comparison
            char csave[60];
            bool lsave[4]  = {false};
            int* isave     = new int[44];
            double* dsave  = new double[29];

            // Function call to LBFGS Fortran subroutine
            setulb_( &n, &mNumCorrections, x, l, u, nbd,
                     &f, g, &mNormDrop, &mGradTolerance,
                     wa, iwa, task, &iprint, csave,
                     lsave, isave, dsave );

            mOptIter = isave[29]; // initiate optimization iteration counter

            while ( (strcmp(task, fg_start) == 0) || (strcmp(task, "NEW_X") == 0) )
            {
                if (strcmp(task, fg_start) == 0)
                {
                    // call to compute objective
                    this->func( mOptIter, x, f );

                    // call to compute gradients
                    this->grad( x, g );
                }

                setulb_( &n, &mNumCorrections, x, l, u, nbd,
                         &f, g, &mNormDrop, &mGradTolerance,
                         wa, iwa, task, &iprint, csave,
                         lsave, isave, dsave );

                std::cout << task << std::endl;

                if ( strcmp(task,"NEW_X") == 0 ) // one iteration of algorithm has concluded
                {
                    mOptIter = isave[29]; // update optimization iteration counter
                }

                // exit loop if maximum iterations have been achieved
                if ( isave[29] == mMaxIt ) break;
            }

            printresult(); // print the result of the optimization algorithm

            aOptProb = mProblem; // update aOptProb

            delete[] nbd;
            delete[] g;
            delete[] wa;
            delete[] iwa;
            delete[] isave;
            delete[] dsave;

            uint tOptIter = gLogger.get_opt_iteration();

            return tOptIter;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_LBFGS::func(int aIter, double* aAdv, double& aObjval )
        {
            // Update the ADV matrix
            Matrix<DDRMat> tADVs(aAdv, mProblem->get_num_advs(), 1);

            // Write restart file
            this->write_advs_to_file(tADVs);

            // Recruit help from other procs and solve for criteria
            this->compute_design_criteria(tADVs);

            // Convert outputs from type MORIS
            aObjval = mProblem->get_objectives()(0);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_LBFGS::grad(double* aAdv, double* aD_Obj )
        {
            auto tD_Obj = mProblem->get_objective_gradients().data();
            std::copy(tD_Obj, tD_Obj + mProblem->get_num_advs(), aD_Obj);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_LBFGS::printresult( )
        {
            std::fprintf( stdout," \nResult of LBFGS\n" );

            // THIS IS A PLACEHOLDER
        }

        //--------------------------------------------------------------------------------------------------------------
    }
}

