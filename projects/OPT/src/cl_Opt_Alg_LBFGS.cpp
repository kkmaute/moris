// Project header files
#include "cl_Opt_Alg_LBFGS.hpp" // OPT/src

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
        OptAlgLBFGS::OptAlgLBFGS( ) :
            OptAlg(),
            mOptIter(0)
        {
            mParameterList.insert( "max_its"  , 100    ); // maximum optimization iterations allowed
            mParameterList.insert( "num_corr" , 5      ); // number of limited memory corrections used in the BFGS update
            mParameterList.insert( "norm_drop", 1.0e+7 ); // LBFGS convergence criteria (this is internally multiplied with machine precision)
            mParameterList.insert( "grad_tol" , 0.0    ); // LBFGS convergence criteria based on projected gradients
        }

        // ---------------------------------------------------------------------

        OptAlgLBFGS::~OptAlgLBFGS()
        {
        }

        //----------------------------------------------------------------------

        void OptAlgLBFGS::solve( OptProb & aOptProb )
        {
            mOptProb = aOptProb;  // set the member variable mOptProb to aOptProb

            OptAlg::initialize(); // initialize the base class member variables

            // extract the underlying types of the algorithm parameters and assign
            // to parameters used by the L-BFGS-B algorithm
            sint tMaxIt  = mParameterList.get< sint >( "max_its" );
            int m        = mParameterList.get< sint >( "num_corr" );
            double factr = mParameterList.get< real >( "norm_drop" );
            double pgtol = mParameterList.get< real >( "grad_tol" );

            int n = mNumAdv; // number of design variables

            // Note that these pointers are deleted by the the Arma and Eigen
            // libraries themselves.
            auto x = mAdvVec.data();
            auto l = mAdvLowVec.data();
            auto u = mAdvUpVec.data();

            // This algorithm parameter is hard-coded to assume that the variables are bounded
            int* nbd = new int[mNumAdv];
            std::fill( nbd, nbd+mNumAdv, 2 );

            // initialize and allocate memory for algorithm inputs/outputs
            double f       = 0;
            double* g      = new double[mNumAdv];
            double* wa     = new double[(2*m + 5)*n + 11*m*m + 8*m];
            int* iwa       = new int[3*mNumAdv];
            int iprint     = -1;      // Prevents the algorithm from printing anything
            char task[61]  = "START"; // starts the algorithm
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

            while ( (strcmp(task,"FG") == 0) || (strcmp(task,"NEW_X") == 0) )
            {
                if (strcmp(task,"FG") == 0)
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

                if ( strcmp(task,"NEW_X") == 0 ) // one iteration of algorithm has concluded
                {
                    mOptIter = isave[29]; // update optimization iteration counter
                }

                // exit loop if maximum iterations have been achieved
                if ( isave[29] == tMaxIt ) break;
            }

            printresult(); // print the result of the optimization algorithm

            aOptProb = mOptProb; // update aOptProb

            delete[] nbd;
            delete[] g;
            delete[] wa;
            delete[] iwa;
            delete[] isave;
            delete[] dsave;
        }

        //----------------------------------------------------------------------

        void OptAlgLBFGS::func( int aIter, double* aAdv, double& aObjval )
        {
            // need to convert to type MORIS to keep coding style consistent
            uint Iter = aIter;

            // Update the ADV matrix
            mAdvVec = Matrix< DDRMat > (aAdv, mNumAdv, 1 );

            // Call to compute objectives and constraints
            OptAlg::func( Iter );

            // Convert outputs from type MORIS
            aObjval = mObjVal;
        }

        //----------------------------------------------------------------------

        void OptAlgLBFGS::grad( double* aAdv, double* aD_Obj )
        {
            // Call to compute derivatives of objectives and constraints
            // w.r.t. advs
            OptAlg::grad( );

            auto tD_Obj = mDObj.data();
            std::copy(tD_Obj, tD_Obj + mNumAdv, aD_Obj);
        }

        //----------------------------------------------------------------------

        void OptAlgLBFGS::printresult( )
        {
            std::fprintf( stdout," \nResult of LBFGS\n" );

            // THIS IS A PLACEHOLDER
        }
    }
}
