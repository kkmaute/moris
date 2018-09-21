// Project header files
#include "cl_Opt_Alg.hpp" // OPT/src

// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        OptAlg::OptAlg() :
                mNumAdv(0),
                mNumCon(0),
                mObjVal(0.0)
        {
        }

        // ---------------------------------------------------------------------

        OptAlg::~OptAlg()
        {
        }

        //----------------------------------------------------------------------

        void OptAlg::initialize()
        {
            mNumAdv   = mOptProb.get_num_adv();   // get the number of advs
            mNumCon   = mOptProb.get_num_con();   // get the number of constraints
            mNumEqCon = mOptProb.get_num_eqcon(); // get the number of equality constraints

            mAdvVec    = mOptProb.get_adv_vec();     // get the advs
            mAdvUpVec  = mOptProb.get_adv_vec_up();  // get the upper bounds of advs
            mAdvLowVec = mOptProb.get_adv_vec_low(); // get the lower bounds of advs

            // set the sizes of the member variable matrices
            mConVal.set_size( mNumCon, 1, 0.0 );       // set the size of the vector of constraints
            mConType.set_size( mNumCon, 1, 1 );        // set the size of the vector of type of constraints
            mDObj.set_size( mNumAdv , 1, 0.0 );        // set the size of the vector of gradient of the objective
            mDCon.set_size( mNumCon , mNumAdv, 0.0 );  // set the size of the matrix of gradient of constraints
            mActive.set_size( mNumCon, 1, 0.0 );       // set the size of the vector containing active constraints flag
        }

        //----------------------------------------------------------------------

        void OptAlg::func( uint aIter )
        {
            // call to OptProb to compute objective and constraints based on
            // the updated design variables
            mOptProb.compute_obj_con( aIter, mAdvVec, mObjVal, mConVal );

            // get the type of constraint flags
            if ( aIter == 0 )
                mConType = mOptProb.get_typecon_vec();
        }

        //----------------------------------------------------------------------

        void OptAlg::grad()
        {
            // call to OptProb to compute gradients of  objective and
            // constraints based on the updated design variables
            mOptProb.grad_obj_con( mActive, mDObj, mDCon );
        }

        //----------------------------------------------------------------------

        void OptAlg::order_con()
        {
            // order constraints such that  1. equality constraints   (typ=0)
            //                              2. inequality constraints (typ=1)

            uint tEqc = 0; // initialize equality constraints counter
            uint tIeq = 0; // initialize inequality constraints counter

            Matrix< DDRMat >  tConVal( mNumCon, 1, 0.0 );

            for ( uint i=0; i<mNumCon; ++i )
            {
                if ( mConType(i,0) == 0 )
                {
                    tConVal( tEqc, 0 ) = mConVal(i,0);
                    tEqc++;
                }
                else
                {
                    tConVal( tIeq+mNumEqCon, 0 ) = mConVal(i,0);
                    tIeq++;
                }
            }

            for ( uint i=0; i<mNumCon; ++i )
            {
                mConVal(i,0) = tConVal(i,0) ;
            }
        }

        //----------------------------------------------------------------------

        void OptAlg::order_grad_con()
        {
            // order constraints such that  1. equality constraints   (typ=0)
            //                              2. inequality constraints (typ=1)

            uint tEqc = 0; // initialize equality constraints counter
            uint tIeq = 0; // initialize inequality constraints counter

            Matrix< DDRMat >  tDCon( mNumAdv, 1, 0.0 );

            for ( uint i=0; i<mNumCon; ++i )
            {
                if ( mConType(i,0) == 0 )
                {
                    tDCon.get_row( tEqc ) = mDCon.get_row(i);
                    tEqc++;
                }
                else
                {
                    tDCon.get_row( tIeq+mNumEqCon ) = mDCon.get_row(i);
                    tIeq++;
                }
            }

            for ( uint i=0; i<mNumCon; ++i )
            {
                mDCon.get_row(i) = tDCon.get_row(i);
            }
        }
    }
}
