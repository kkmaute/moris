// Project header files
#include "cl_OPT_Algorithm_LBFGS.hpp" // OPT/src

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
            double* f, double* g, double* factr, double* pgtol,
            double* wa, int* iwa, char* task, int* iprint, char* csave,
            bool* lsave, int* isave, double* dsave);
}

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Algorithm_LBFGS::Algorithm_LBFGS( ) :
                Algorithm(),
                mOptIter(0)
        {
            mParameterList.insert( "max_its"  , 100    ); // maximum optimization iterations allowed
            mParameterList.insert( "num_corr" , 5      ); // number of limited memory corrections used in the BFGS update
            mParameterList.insert( "norm_drop", 1.0e+7 ); // LBFGS convergence criteria (this is internally multiplied with machine precision)
            mParameterList.insert( "grad_tol" , 0.0    ); // LBFGS convergence criteria based on projected gradients
        }

        // ---------------------------------------------------------------------

        Algorithm_LBFGS::~Algorithm_LBFGS()
        {
        }

        //----------------------------------------------------------------------

        void Algorithm_LBFGS::solve(Problem* aOptProb )
        {
            mProblem = aOptProb;  // set the member variable mProblem to aOptProb

            Algorithm::initialize(); // initialize the base class member variables

            // extract the underlying types of the algorithm parameters and assign
            // to parameters used by the L-BFGS-B algorithm
            sint tMaxIt  = mParameterList.get< sint >( "max_its" );
            int m        = mParameterList.get< sint >( "num_corr" );
            double factr = mParameterList.get< real >( "norm_drop" );
            double pgtol = mParameterList.get< real >( "grad_tol" );

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
            double* wa     = new double[(2*m + 5)*n + 11*m*m + 8*m];
            int* iwa       = new int[3 * n];
            int iprint     = -1;      // Prevents the algorithm from printing anything
            char task[61]  = "START"; // starts the algorithm
            char fg_start[61] = "FG_START                                                    "; // needed for comparison
            char csave[60];
            bool lsave[4]  = {false};
            int* isave     = new int[44];
            double* dsave  = new double[29];

            // Function call to LBFGS Fortran subroutine
            setulb_( &n, &m, x, l, u, nbd,
                     &f, g, &factr, &pgtol,
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

                setulb_( &n, &m, x, l, u, nbd,
                         &f, g, &factr, &pgtol,
                         wa, iwa, task, &iprint, csave,
                         lsave, isave, dsave );

                std::cout << task << std::endl;

                if ( strcmp(task,"NEW_X") == 0 ) // one iteration of algorithm has concluded
                {
                    mOptIter = isave[29]; // update optimization iteration counter
                }

                // exit loop if maximum iterations have been achieved
                if ( isave[29] == tMaxIt ) break;
            }

            printresult(); // print the result of the optimization algorithm

            aOptProb = mProblem; // update aOptProb

            delete[] nbd;
            delete[] g;
            delete[] wa;
            delete[] iwa;
            delete[] isave;
            delete[] dsave;
        }

        //----------------------------------------------------------------------

        void Algorithm_LBFGS::func(int aIter, double* aAdv, double& aObjval )
        {
            // Update the ADV matrix
            Matrix<DDRMat> tADVs(aAdv, mProblem->get_num_advs(), 1);
            mProblem->set_advs(tADVs);

            // Call to update objectives and constraints
            mProblem->mUpdateObjectives = true;

            // Convert outputs from type MORIS
            aObjval = mProblem->get_objectives()(0);
        }

        //----------------------------------------------------------------------

        void Algorithm_LBFGS::grad(double* aAdv, double* aD_Obj )
        {
            // Call to update derivatives of objectives
            mProblem->mUpdateObjectiveGradient = true;

            auto tD_Obj = mProblem->get_objective_gradients().data();
            std::copy(tD_Obj, tD_Obj + mProblem->get_num_advs(), aD_Obj);
        }

        //----------------------------------------------------------------------

        void Algorithm_LBFGS::printresult( )
        {
            std::fprintf( stdout," \nResult of LBFGS\n" );

            // THIS IS A PLACEHOLDER
        }
    }
}
