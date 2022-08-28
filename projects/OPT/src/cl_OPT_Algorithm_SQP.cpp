/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_SQP.cpp
 *
 */

#include "cl_OPT_Algorithm_SQP.hpp"
#include "cl_Communication_Tools.hpp"

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
// SNOPT (a Fortran implementation of SQP) function declarations
void sninit_(
        int* tPrint, int* tSumm, char* cw, int* lencw,
        int* iw, int* leniw, double* rw, int* lenrw);

void snset_(
        char* buffer, int* tPrint, int* tSumm, int* tExit,
        char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
        short buflen);

void snseti_(
        char* buffer, int* ivalue, int* tPrint, int* tSumm, int* tExit,
        char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
        short buflen);

void sngeti_(
        char* buffer, int* ivalue, int* tExit,
        char* cw, int* lencw, int* iw, int* leniw, double* rw, int* lenrw,
        short buflen);

void snsetr_(
        char* buffer, double* rvalue, int* tPrint, int* tSumm, int* tExit,
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
        int* tInform, int* mincw, int* miniw, int* minrw,
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

// ---------------------------------------------------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        Algorithm_SQP::Algorithm_SQP(ParameterList aParameterList)
        {
            // Initialize
            int tPrint = 0;
            int tSumm  = 0;
            int tExit  = 0;

            // Remove algorithm name from the parameter list
            aParameterList.erase("algorithm");

            // loop over all entries in the parameter list
            for ( auto it = aParameterList.begin(); it != aParameterList.end(); ++it )
            {
                // skip over non-SQP specific parameters
                if ( it->first == "restart_index") continue;

                sint tParamType = aParameterList.which( it->first ); // get the underlying parameter type

                // call Fortran subroutine based on parameter type
                switch ( tParamType )
                {
                    case 1: // set integer parameters
                    {
                        char* paramname = (char*)it->first.c_str();
                        int tParamVal   = aParameterList.get< sint >( paramname );
                        short paramlen  = strlen( paramname );

                        snseti_( paramname, &tParamVal, &tPrint, &tSumm, &tExit,
                                mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW, paramlen );

                        break;
                    }

                    case 2: // set double parameters
                    {
                        char* paramname  = (char*)it->first.c_str();
                        double tParamVal = aParameterList.get< real >( paramname );
                        short paramlen   = strlen( paramname );

                        snsetr_( paramname, &tParamVal, &tPrint, &tSumm, &tExit,
                                mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW, paramlen );
                        break;
                    }

                    case 4: // set string parameters to char
                    {
                        char* paramname = (char*)it->first.c_str();
                        char* tParamVal = (char*)( aParameterList.get< std::string >( paramname ) ).c_str();
                        short paramlen  = strlen( tParamVal ) ;

                        snset_( tParamVal, &tPrint, &tSumm, &tExit,
                                mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW, paramlen );
                        break;
                    }

                    default:
                        MORIS_LOG_ERROR ( "No matching function call for underlying type.");
                        assert::error( "In cl_Algorithm_SQP.cpp" );
                }

                if( tExit != 0 )
                {
                    MORIS_LOG_ERROR ( "When calling SNOPT Fortran subroutine, unable to set parameter :  %-5i", it->first.c_str());
                    assert::error( "In cl_Algorithm_SQP.cpp" );
                }
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        Algorithm_SQP::~Algorithm_SQP()
        {
        }

        //--------------------------------------------------------------------------------------------------------------

        uint Algorithm_SQP::solve( uint aCurrentOptAlgInd, std::shared_ptr<Problem> aOptProb )
        {
            // Trace optimization
            Tracer tTracer( "OptimizationAlgorithm", "SQP", "Solve" );

            mCurrentOptAlgInd = aCurrentOptAlgInd;  // set index of current optimization algorithm
            mProblem          = aOptProb;           // set the member variable mProblem to aOptProb

            // Set optimization iteration index for restart
            if ( mRestartIndex > 0)
            {
                gLogger.set_opt_iteration( mRestartIndex );
            }

            // Solve optimization problem
            if (par_rank() == 0)
            {
                // Run sqp algorithm
                this->sqp_solve();

                // update aOptProb FIXME: Why is this needed
                aOptProb = mProblem;

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

            return tOptIter;
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_SQP::sqp_solve()
        {
            int tInform = 0;
            int tPrint = 0;
            int tSumm  = 0;

            // initialize
            sninit_(&tPrint, &tSumm, mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW);

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
            //set_params(cw, mLenCW, mIW, mLenIW, mRW, mLenRW);

            // get parameter values that have been changed internally
            {
                char* paramname = (char*)"Print file";

                short paramlen  = strlen(paramname);
                sngeti_(paramname, &tPrint, &tInform,
                        mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW, paramlen);
            }

            {
                char* paramname = (char*)"Summary file";
                short paramlen  = strlen(paramname);
                sngeti_(paramname, &tSumm, &tInform,
                        mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW, paramlen);
            }
            // estimate workspace sizes
            snmema_(&nF, &n, &nxname, &nFname, &neA, &neG,
                    &mLenCW, &mLenIW, &mLenRW, mCW, &mMinWLen, mIW, &mMinWLen, mRW, &mMinWLen);

            mCW = static_cast<char*>  (realloc(mCW, mLenCW*8*sizeof(char)));
            mIW = static_cast<int*>   (realloc(mIW, mLenIW*sizeof(int)));
            mRW = static_cast<double*>(realloc(mRW, mLenRW*sizeof(double)));

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
            double rmem = sizeof(double)*(mLenRW+2*n+4*nF) + sizeof(int)*(nF+n+mLenIW+lenG) + sizeof(char)*(mLenCW*8);

            MORIS_LOG_INFO("Allocating memory for SNOPT: %8.2f Mb.\n",rmem/1024.0/1024.0);

            // initialize with resized work arrays
            sninit_(&tPrint, &tSumm, mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW);

            // set parameters
            //set_params(cw, mLenCW, mIW, mLenIW, mRW, mLenRW);

            // solve problem
            snopta_(&Start,
                    &nF, &n, &nxname, &nFname,
                    &mObjAdd, &mObjRow, mProb,
                    reinterpret_cast<void (*)()>(sqp_user_function),
                    reinterpret_cast<void (*)()>(snlog_),
                    &iAfun, &jAvar, &lenA, &neA, &A,
                    iGfun, jGvar, &lenG, &neG,
                    xlow, xupp, xnames,
                    Flow, Fupp, Fnames,
                    x, xstate, xmul,
                    F, Fstate, Fmul,
                    &tInform, &mincw, &miniw, &minrw,
                    &nS, &nInf, &sInf,
                    cu, &lencu, iu, &leniu, ru, &lenru,
                    mCW, &mLenCW, mIW, &mLenIW, mRW, &mLenRW,
                    strlen(mProb));

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
            free(mCW);
            free(mIW);
            free(mRW);
        }

        //--------------------------------------------------------------------------------------------------------------

        void Algorithm_SQP::func_grad(int n, double* x, int needG )
        {
            // get ADVs from problem
            auto tAdvVec = mProblem->get_advs().data();

            // check whether criteria need to be evaluated
            if( mOptIter == 0 or !std::equal(x, x+n, tAdvVec) )
            {
                // update the vector of design variables
                Matrix<DDRMat> tADVs(x, mProblem->get_num_advs(), 1);

                // Write restart file
                this->write_advs_to_file(tADVs);

                // Compute design criteria
                this->compute_design_criteria(tADVs);

                // Compute design criteria gradients
                this->compute_design_criteria_gradients(tADVs);

                // set update for objectives and constraints
                if (needG)
                {
                    // FIXME figure out a way to move dcriteria_dadv_solve here while still having parallel work
                }

                ++mOptIter;
            }
        }

        //--------------------------------------------------------------------------------------------------------------

        void sqp_user_function(
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
                moris::Matrix< moris::DDRMat > tGradObj = tOptAlgSQP->mProblem->get_objective_gradients();
                moris::Matrix< moris::DDRMat > tGradCon = tOptAlgSQP->mProblem->get_constraint_gradients();

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

        //--------------------------------------------------------------------------------------------------------------

    }
}

