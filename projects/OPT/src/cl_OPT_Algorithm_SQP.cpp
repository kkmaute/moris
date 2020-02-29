// Project header files
#include "cl_OPT_Algorithm_SQP.hpp" // OPT/src
extern moris::Logger gLogger;

#ifdef FORT_NO_
#define _FORTRAN(a) a
#elif  F77ADD_
#define _FORTRAN(a) a ## _
#elif F77ADD__
#define _FORTRAN(a) a ## __
#endif

extern "C"
{
    // SNOPT (a Fortran implementation of SQP) function declarations
    void sninit_(
            int* iPrint, int* iSumm, char* cw, int* lencw,
            int* iw, int* leniw, double* rw, int* lenrw);

    void snset_(
            char* buffer, int* iPrint, int* iSumm, int* iExit,
            char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
            short buflen);

    void snseti_(
            char* buffer, int* ivalue, int* iPrint, int* iSumm, int* iExit,
            char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
            short buflen);

    void sngeti_(
            char* buffer, int* ivalue, int* iExit,
            char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
            short buflen);

    void snsetr_(
            char* buffer, double* rvalue, int* iPrint, int* iSumm, int* iExit,
            char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
            short buflen);

    void snopta_(
            int* Start, int* nF, int* n, int* nxname, int* nFname,
            double* mObjAdd, int* mObjRow, char* mProb, void (*usrfun)(), void (*snLog)(),
            int* iAfun, int* jAvar, int* lenA, int* neA, double* A,
            int* iGfun, int* jGvar, int* lenG, int* neG,
            double* xlow, double* xupp, char* xnames,
            double* Flow, double* Fupp, char* Fnames,
            double* x, int* xstate, double* xmul,
            double* F, int* Fstate, double* Fmul,
            int* inform, int* mincw, int* miniw, int* minrw,
            int* nS, int* nInf, double* sInf,
            char* cu, int* lencu, int* iu, int* leniu,
            double* ru, int* lenru,
            char* cw, int* lencw, int* iw, int* leniw,
            double* rw, int* lenrw, short mProblen);

    void snmema_(
            int* nF, int* n, int* nxname, int* nFname, int* neA, int* neG,
            int* mincw, int* miniw, int* minrw, char* cw, int* lencw, int* iw,
            int* leniw, double* rw, int* lenrw);

    void snlog_(
            int* iAbort, int* info, int* Htype, int* KTcond, int* MjrPrt, int* minimz,
            int* n, int* nb, int* nnCon0, int* nS, int* itn, int* nMajor, int* nMinor, int* nSwap,
            double* condHz, int* iObj, double* sclObj, double* mObjAdd, double* fMrt, double* PenNrm, double* step,
            double* prInf, double* duInf, double* vimax, double* virel, int* hs,
            int* ne, int* nlocJ, int* locJ, int* indJ, double* Jcol,
            double* Ascale, double* bl, double* bu, double* Fcon, double* Lmul, double* x,
            char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru,
            char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw);
}

void OptAlgSQP_usrfun(
        int* Status, int* n, double* x,
        int* needf, int* nF, double* f,
        int* needG, int* lenG, double* G,
        char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru);

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Algorithm_SQP::Algorithm_SQP( ) :
                Algorithm(),
                mOptIter(0)
        {
            // assign algorithm specific parameters - these values have been
            // carried forward from the old code
            mMinWLen = 500;
            mObjAdd  = 0.0;
            mObjRow  = 1;
            mProb    = (char*)"fem ";

            //******************************************************************
            // assign default parameter values - PLEASE REFER to documentation
            // at '.../MTPLS/Snopt/doc/snopt.pdf' for CLEAR UNDERSTANDING of the
            // parameters encountered in this algorithm.
            //******************************************************************

            // Printing
            mParameterList.insert( "Major print level", 1 );   // Controls the amount of output to print each major iteration.
            mParameterList.insert( "Minor print level", 1 );   // Controls the amount of output during the QP subproblems.
            mParameterList.insert( "Print file"       , 0 );   // Change to >0 to print output to file.
            mParameterList.insert( "Summary file"     , 0 );   // Change to >0 to print summary to file.
            mParameterList.insert( "Print frequency"  , 100 ); // Every nth minor iteration we print output to file.
            mParameterList.insert( "Log frequency"    , 100 ); // Related to print frequency.
            mParameterList.insert( "Summary frequency", 100 ); // Every nth minor iteration we print output to file.
            mParameterList.insert( "Timing level"     , 3 );   // prints CPU times

            // SQP method
            mParameterList.insert( "Major iterations limit", 1000 );                    // number of allowed major iterations
            mParameterList.insert( "Minor iterations limit", 500 );                     // number of allowed minor iterations
            mParameterList.insert( "Major step limit"      , 2.0 );                     // limits the change in variable during linesearch
            mParameterList.insert( "Superbasics limit"     , 500 );                     // places a limit on the storage of superbasic variables
            mParameterList.insert( "New superbasics limit" , 99 );                      // controls early termination of QPs
            mParameterList.insert( "linesearch_type"       , "Derivative linesearch" ); // other options are:
                                                                                        // "Nonderivative linesearch"
            mParameterList.insert( "Linesearch tolerance"       , 0.9 );                // controls accuracy of linesearch
            mParameterList.insert( "Function precision"         , 3e-13 );              // relative accuracy with which nonlinear functions are computed
            mParameterList.insert( "Difference interval"        , 5.5e-7 );             // determines accuracy of gradients using forward differencing
            mParameterList.insert( "Central difference interval", 6.7e-5 );             // determines accuracy of gradients using central differencing
            mParameterList.insert( "Proximal point method"      , 1 );                  // satisfies linear constraints near initial guess
            mParameterList.insert( "Violation limit"            , 10.0 );               // limit on maximum constraint violation after linesearch
            mParameterList.insert( "Unbounded step size"        , 1.0e18 );             // determines unboundedness of linesearch step size
            mParameterList.insert( "Unbounded objective value"  , 1.0e15 );             // determines unboundedness of objective
            mParameterList.insert( "Infinite bound size", 1.0e+20 );                    // any upper bound greater than this value is regarded as infinity

            // QP subproblems
            mParameterList.insert( "Iterations limit", 10000 );   // sum of allowed minor iterations over all major iterations
            mParameterList.insert( "Elastic weight"  , 2.0e+4 );  // weighting of infeasibilities in the objective of the QP subproblem
            mParameterList.insert( "Partial price"   , 1 );       // reduces the work required for each "pricing" operation
            mParameterList.insert( "Pivot tolerance" , 3.7e-11 ); // guards the basis matrix from becoming singular

            // Hessian approximation
            mParameterList.insert( "hessian_type"     , "Hessian Full memory" ); // Method for storing and updating the Hessian.
                                                                                 // Set to "Hessian Limited memory" for variables > 75.
            mParameterList.insert( "Hessian frequency", 999999 );                // for full memory Hessian
            mParameterList.insert( "Hessian updates"  , 20 );                    // for limited memory Hessian

            // Frequencies
            mParameterList.insert( "Expand frequency"       , 10000 ); // for anti-cycling procedure
            mParameterList.insert( "Factorization frequency", 50 );    // for basis updates

            // LU options
            mParameterList.insert( "LU factor tolerance"     , 10.0 );                  // limits size of multipliers in L
            mParameterList.insert( "LU update tolerance"     , 10.0 );                  // limits size of multipliers in L during updates
            mParameterList.insert( "LU density tolerance"    , 0.6  );                  // handles sparsity of LU factorization
            mParameterList.insert( "LU singularity tolerance", 2e-6 );                  // handles guard against singularity during factorization
            mParameterList.insert( "lu_pivoting_type"        , "LU Partial Pivoting" ); // Related to LU factorization. Other options are
                                                                                        // "LU Rook Pivoting" - more costly and stable
                                                                                        // "LU Complete Pivoting" - more costly and stable

            // Convergence Tolerances
            mParameterList.insert( "Major optimality tolerance" , 1e-6 ); // target accuracy of the dual variable
            mParameterList.insert( "Minor optimality tolerance" , 5e-7 ); // Also related to target accuracy of the dual variable
            mParameterList.insert( "Major feasibility tolerance", 1e-6 ); // target nonlinear constraint violation
            mParameterList.insert( "Feasibility tolerance"      , 1e-6 ); // Related to minor feasibility tolerance
            mParameterList.insert( "Minor feasibility tolerance", 1e-6 ); // for satisfying the QP bounds

            // Derivative checking
            mParameterList.insert( "Verify level", 0 ); // Finite difference check on derivatives computed by user-provided routines

            // Scaling
            mParameterList.insert( "Scale option"   , 1 );   // flag for scaling of constraints and variables
            mParameterList.insert( "Scale tolerance", 0.9 ); // affects effectiveness with which constraint matrix is scaled
        }

        //----------------------------------------------------------------------

        Algorithm_SQP::~Algorithm_SQP()
        {
        }

        //----------------------------------------------------------------------

        void Algorithm_SQP::set_params(char* cw, int lencw, int* iw, int leniw, double* rw, int lenrw)
        {
            int iPrint = 0;
            int iSumm  = 0;
            int iExit  = 0;

            // loop over all entries in the parameter list
            for ( auto it = mParameterList.begin(); it != mParameterList.end(); ++it )
            {
                sint tParamType = mParameterList.which( it->first ); // get the underlying parameter type

                // call Fortran rubroutine based on parameter type
                switch ( tParamType )
                {
                case 1: // set integer parameters
                {
                    char* paramname = (char*)it->first.c_str();
                    int tParamVal   = mParameterList.get< sint >( paramname );
                    short paramlen  = strlen( paramname );

                    snseti_( paramname, &tParamVal, &iPrint, &iSumm, &iExit,
                             cw, &lencw, iw, &leniw, rw, &lenrw, paramlen );

                    break;
                }

                case 2: // set double parameters
                {
                    char* paramname  = (char*)it->first.c_str();
                    double tParamVal = mParameterList.get< real >( paramname );
                    short paramlen   = strlen( paramname );

                    snsetr_( paramname, &tParamVal, &iPrint, &iSumm, &iExit,
                             cw, &lencw, iw, &leniw, rw, &lenrw, paramlen );
                    break;
                }

                case 4: // set string parameters to char
                {
                    char* paramname = (char*)it->first.c_str();
                    char* tParamVal = (char*)( mParameterList.get< std::string >( paramname ) ).c_str();
                    short paramlen  = strlen( tParamVal ) ;

                    snset_( tParamVal, &iPrint, &iSumm, &iExit,
                            cw, &lencw, iw, &leniw, rw, &lenrw, paramlen );
                    break;
                }

                default:
                    MORIS_LOG_ERROR ( "No matching function call for underlying type.");
                    assert::error( "In cl_Algorithm_SQP.cpp" );
                }

                if( iExit != 0 )
                {
                    MORIS_LOG_ERROR ( "When calling SNOPT Fortran subroutine, unable to set parameter :  %-5i", it->first.c_str());
                    assert::error( "In cl_Algorithm_SQP.cpp" );
                }
            }
        }

        //----------------------------------------------------------------------

        void Algorithm_SQP::solve(Problem* aOptProb )
        {
            mProblem = aOptProb; // set the member variable mProblem to aOptProb

            Algorithm::initialize(); // initialize the base class member variables

            int inform = 0;
            int iPrint = 0;
            int iSumm  = 0;

            int  lencw  = mMinWLen;
            char* cw    = static_cast<char*>  (malloc(lencw*8*sizeof(char)));

            int  leniw  = mMinWLen;
            int* iw     = static_cast<int*>   (malloc(leniw*sizeof(int)));

            int  lenrw  = mMinWLen;
            double* rw  = static_cast<double*>(malloc(lenrw*sizeof(double)));

            // initialize
            sninit_(&iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw);

            const int nvar  = mProblem->get_num_advs();    // number of abstract variables
            const int m_con = mProblem->get_num_constraints();    // total number of constraints
            const int m_eq  = mProblem->get_num_equality_constraints();  // number of equality constraints

            int nF     = m_con + 1;
            int n      = nvar;
            int nxname = 1;
            int nFname = 1;
            int neA    = 0;
            int neG    = n*nF;
            int lenG   = neG;

            int* iGfun = static_cast<int*>(malloc(lenG*sizeof(int)));
            int* jGvar = static_cast<int*>(malloc(lenG*sizeof(int)));

            for(int i = 0, cur_nz = 0; i < n; ++i)
            {
                for(int j = 0; j < nF; ++j, ++cur_nz)
                {
                    iGfun[cur_nz] = j+1;
                    jGvar[cur_nz] = i+1;
                }
            }

            int lenA = 1;
            int iAfun, jAvar;
            double A;

            // extract design variables including their lower and upper bounds
            auto x    = mProblem->get_advs().data();
            auto xlow = mProblem->get_lower_bounds().data();
            auto xupp = mProblem->get_upper_bounds().data();

            mProblem->get_advs().fill(MAXDOUBLE); // The overall logic needs this

            char*   Fnames = 0;
            char*   xnames = 0;

            double* Flow = static_cast<double*>(malloc(nF*sizeof(double)));
            double* Fupp = static_cast<double*>(malloc(nF*sizeof(double)));

            std::fill(Flow       , Flow+1     , -MAXDOUBLE);
            std::fill(Flow+1     , Flow+nF    , 0.0);
            std::fill(Fupp       , Fupp+1     , MAXDOUBLE);
            std::fill(Fupp+1     , Fupp+1+m_eq, 0.0);
            std::fill(Fupp+1+m_eq, Fupp+nF    , MAXDOUBLE);

            double* xmul = static_cast<double*>(malloc(n*sizeof(double)));
            double* F    = static_cast<double*>(malloc(nF*sizeof(double)));
            double* Fmul = static_cast<double*>(malloc(nF*sizeof(double)));

            int*   Fstate = static_cast<int*>(malloc(nF*sizeof(int)));
            int*   xstate = static_cast<int*>(malloc(n*sizeof(int)));

            std::fill(xstate, xstate+n, 0);

            // set parameters
            set_params(cw, lencw, iw, leniw, rw, lenrw);

            // get parameter values that have been changed itnernally
            {
                char* paramname = (char*)"Print file";
                short paramlen  = strlen(paramname);
                sngeti_(paramname, &iPrint, &inform,
                        cw, &lencw, iw, &leniw, rw, &lenrw, paramlen);
            }
            {

                char* paramname = (char*)"Summary file";
                short paramlen  = strlen(paramname);
                sngeti_(paramname, &iSumm, &inform,
                        cw, &lencw, iw, &leniw, rw, &lenrw, paramlen);
            }

            // estimate workspace sizes
            snmema_(&nF, &n, &nxname, &nFname, &neA, &neG,
                    &lencw, &leniw, &lenrw, cw, &mMinWLen, iw, &mMinWLen, rw, &mMinWLen);

            cw = static_cast<char*>  (realloc(cw, lencw*8*sizeof(char)));
            iw = static_cast<int*>   (realloc(iw, leniw*sizeof(int)));
            rw = static_cast<double*>(realloc(rw, lenrw*sizeof(double)));

            int mincw, miniw, minrw, nS, nInf;
            double sInf;

            int lencu  = 0;
            int leniu  = 0;
            int lenru  = 0;
            int Start  = 0;
            int* iu    = 0;
            double* ru = 0;

            char* cu = reinterpret_cast<char*>(this);

            // calculate required memory
            double rmem = sizeof(double)*(lenrw+2*n+4*nF) + sizeof(int)*(nF+n+leniw+lenG) + sizeof(char)*(lencw*8);

            if( iPrint > 0 )
                fprintf(stderr," ... Allocating memory for SNOPT: %8.2f Mb.\n",rmem/1024.0/1024.0);

            // initialize with resized work arrays
            sninit_(&iPrint, &iSumm, cw, &lencw, iw, &leniw, rw, &lenrw);

            // set parameters
            set_params(cw, lencw, iw, leniw, rw, lenrw);

            // solve mProblem
            snopta_(&Start,
                    &nF, &n, &nxname, &nFname,
                    &mObjAdd, &mObjRow, mProb,
                    reinterpret_cast<void (*)()>(OptAlgSQP_usrfun),
                    reinterpret_cast<void (*)()>(snlog_),
                    &iAfun, &jAvar, &lenA, &neA, &A,
                    iGfun, jGvar, &lenG, &neG,
                    xlow, xupp, xnames,
                    Flow, Fupp, Fnames,
                    x, xstate, xmul,
                    F, Fstate, Fmul,
                    &inform, &mincw, &miniw, &minrw,
                    &nS, &nInf, &sInf,
                    cu, &lencu, iu, &leniu, ru, &lenru,
                    cw, &lencw, iw, &leniw, rw, &lenrw,
                    strlen(mProb));

            aOptProb = mProblem; // update aOptProb

            // delete temporary arrays
            free(F);
            free(Fstate);
            free(Fmul);
            free(xstate);
            free(xmul);
            free(Flow);
            free(Fupp);
            free(iGfun);
            free(jGvar);
            free(cw);
            free(iw);
            free(rw);
        }

        //----------------------------------------------------------------------

        void Algorithm_SQP::func_grad(int n, double* x, int needG )
        {
            auto tAdvVec = mProblem->get_advs().data();

            if( !std::equal(x, x+n, tAdvVec) )
            {
                // update the vector of design variables
                Matrix<DDRMat> tADVs(x, mProblem->get_num_advs(), 1);
                mProblem->set_advs(tADVs);

                // set update for objectives and constraints
                mProblem->mUpdateObjectives = true;

                if(needG)
                {
                    mProblem->mUpdateConstraints = true;
                }

                ++mOptIter;
            }
        }
    }
}

//----------------------------------------------------------------------

void OptAlgSQP_usrfun(
        int* Status, int* n, double* x,
        int* needf, int* nF, double* f,
        int* needG, int* lenG, double* G,
        char* cu, int* lencu, int* iu, int* leniu, double* ru, int* lenru )
{
    moris::opt::Algorithm_SQP* tOptAlgSQP = reinterpret_cast< moris::opt::Algorithm_SQP* >(cu);
    tOptAlgSQP->func_grad(*n, x, *needG);
    *Status = 0;

    if(*needf)
    {
        f[0] = tOptAlgSQP->mProblem->get_objectives()(0); // obtain the objective

        // obtain the vector of constraints
        moris::Matrix< moris::DDRMat > tConVal = tOptAlgSQP->mProblem->get_constraints();

        // Convert constraints to a pointer from moris::mat
        for ( moris::uint i=0; i < tConVal.n_rows(); ++i )
        {
            // tConVal is of type MORIS matrix
            f[1+i] = -tConVal(i,0);  // Need a negative sign such that the definition of constraints is consistent in the input file
        }
    }

    if(*needG)
    {
        // get the gradient of objective and constraints
        moris::Matrix< moris::DDRMat > tGradObj = tOptAlgSQP->mProblem->get_objective_gradient();
        moris::Matrix< moris::DDRMat > tGradCon = tOptAlgSQP->mProblem->get_constraint_gradient();

        for( int i = 0, cur_nz = 0; i < *n; ++i )
        {
            G[cur_nz++] = tGradObj(i);

            for( int j = 1; j < *nF; ++j, ++cur_nz )
            {
                G[cur_nz] = -tGradCon(j-1, i); // Need a negative sign such that the definition of constraints is consistent in the input file
            }
        }
    }
}
