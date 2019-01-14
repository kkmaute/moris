// Project header files
#include "cl_Opt_Prob.hpp" // OPT/src
#include "op_plus.hpp"
#include "Opt_Input.hpp"
extern moris::Logger gLogger;
// -----------------------------------------------------------------------------

namespace moris
{
    namespace opt
    {
        OptProb::OptProb(
                uint aNumMasterAdvs,
                uint aNumCriteria,
                uint aNumConstraints,
                uint aNumEqConstraints ) :
            mNumMasterAdvs( aNumMasterAdvs ),
            mNumCriteria( aNumCriteria ),
            mNumConstraints( aNumConstraints ),
            mNumEqConstraints( aNumEqConstraints ),
            mObj( 0.0 )
        {
            // set the size of all the member variable matrices and initialize
            // with zeros
            mAbsDesVarVec.set_size( mNumMasterAdvs, 1, 0.0 );;         // Abstract Design Variable Vector
            mAbsDesVarVecUp.set_size( mNumMasterAdvs, 1, 0.0 );;       // Upper bounds on Abstract Design Variable Vector
            mAbsDesVarVecLow.set_size( mNumMasterAdvs, 1, 0.0 );;      // Lower bounds on  Abstract Design Variable Vector
            mCon.set_size( mNumConstraints, 1, 0.0 );                  // matrix of constraints
            mTypeCon.set_size( mNumConstraints, 1, -1 );               // matrix containing flags for type of constraints
            mOptCriteria.set_size( mNumCriteria, 1, 0.0 );             // Set of optimization criteria
            mGradObj.set_size( mNumMasterAdvs, 1, 0.0 );               // derivative of objective w.r.t advs
            mGradCon.set_size( mNumConstraints, mNumMasterAdvs, 0.0 ); // derivative of constraints w.r.t advs
            mGradCrit.set_size( mNumCriteria, mNumMasterAdvs, 0.0 );   // derivative of criteria w.r.t advs
            mActive.set_size( mNumConstraints, 1, 0 );                 // flag for active/inactive constraints

            this->define_advs();         // user defined initial advs
            this->define_opt_criteria(); // user defined optimization criteria
        }

        // ---------------------------------------------------------------------

        OptProb::~OptProb()
        {
        }

        // ---------------------------------------------------------------------

        void OptProb::define_advs()
        {
            // obtain the advs and their upper and lower bound from the input
            // file
            opt_input::define_advs( mAbsDesVarVec, mAbsDesVarVecUp, mAbsDesVarVecLow );
        }

        // ---------------------------------------------------------------------

        void OptProb::define_opt_criteria()
        {
            // obtain the optimization criteria definition from the input file
            opt_input::define_opt_criteria( mOptCriteria );
        }

        // ---------------------------------------------------------------------

        void OptProb::compute_obj_con(
                const uint          aIter,
                const Matrix< DDRMat >  & aAdvVec,
                real              & aObjVal,
                Matrix< DDRMat >        & aConVal )
        {
            mAbsDesVarVec = aAdvVec; // update the adv vector

            // The physics is decided based upon the criteria opted for
            this->compute_criteria();

            // evaluate the objective and constraints based on their definition
            // in the input file
            opt_input::get_obj_con( aIter, mAbsDesVarVec, mOptCriteria, mObj, mCon, mTypeCon );

            aObjVal = mObj; // assign the updated objective value
            aConVal = mCon; // assign the updated constraint values
        }

        // ---------------------------------------------------------------------

        void OptProb::compute_criteria()
        {
            // Call to a more generic compute scalar would come over here
            // The criteria would usually be something more complicated in
            // presence of the FEA module

            if ( mIsOptAnlyt )
            {
                // Get the criteria from the input file
                opt_input::get_criteria( mAbsDesVarVec, mOptCriteria );
            }
            else
            {
                MORIS_LOG_ERROR ( "Only analytical criteria possible.\n");
                assert::error( "In cl_Opt_Prob.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        // Once we have the concept of dofs, this function would have an option
        // to compute du/ds = J(^-1)dR/ds using the adjoint method, direct
        // method or through finite difference

        void OptProb::grad_obj_con(
                const Matrix< DDSMat >  & aActive,
                Matrix< DDRMat >        & aGradObj,
                Matrix< DDRMat >        & aGradCon )
        {
            mActive = aActive; // update the matrix of active constraint flags

            // create matrices for derivative of objective and constraints
            Matrix< DDRMat >  tGradObj_ds( mNumMasterAdvs, 1.0, 0.0 );             // matrix for explicit gradient of objective w.r.t. adv
            Matrix< DDRMat >  tGradCon_ds( mNumConstraints, mNumMasterAdvs, 0.0 );  // matrix for explicit gradient of constraints w.r.t. adv
            Matrix< DDRMat >  tGradObj_dCrit( mNumCriteria, 1.0, 0.0 );          // matrix for explicit gradient of objective w.r.t. criteria
            Matrix< DDRMat >  tGradCon_dCrit( mNumConstraints, mNumCriteria, 0.0 ); // matrix for explicit gradient of constraints w.r.t. criteria

            // Compute explicit derivative of objective and constraints w.r.t.
            // criteria
            this->explicit_dobjcon_dcrit( tGradObj_dCrit, tGradCon_dCrit );

            // Compute derivative of criteria w.r.t. adv
            this->dcriteria_ds();

            // Compute explicit derivative of objective and constraints w.r.t. adv
            this->explicit_dobjcon_ds( tGradObj_ds, tGradCon_ds );

            // Assemble the overall objective and constraint gradients
            this->assemble_grad( tGradObj_dCrit, tGradCon_dCrit, tGradObj_ds, tGradCon_ds );

            aGradObj = mGradObj; // assign the objective gradient value
            aGradCon = mGradCon; // assign the constraint gradient values
        }

        // ---------------------------------------------------------------------

        void OptProb::explicit_dobjcon_dcrit(
                Matrix< DDRMat >  & aGradObj_dCrit,
                Matrix< DDRMat >  & aGradCon_dCrit )
        {
            // Get the explicit gradients of objective and constraints w.r.t
            // optimization criteria from the input file
            opt_input::get_dobjcon_dcrit( mAbsDesVarVec, mOptCriteria, aGradObj_dCrit, aGradCon_dCrit );
        }

        // ---------------------------------------------------------------------

        void OptProb::dcriteria_ds()
        {
            // Call to a more generic routine would come over here
            // The criteria would usually be something more complicated in
            // presence of the FEA module

            if ( mIsOptAnlyt )
            {
                // Get the derivative of the criteria from the input file
                opt_input::get_criteria( mAbsDesVarVec, mGradCrit );
            }
            else
            {
                MORIS_LOG_ERROR ( "Only analytical criteria possible.\n");
                assert::error( "In cl_Opt_Prob.cpp" );
            }
        }

        // ---------------------------------------------------------------------

        void OptProb::explicit_dobjcon_ds(
                Matrix< DDRMat >  & aGradObj_ds,
                Matrix< DDRMat >  & aGradCon_ds)
        {
            // Get the explicit gradients of objective and constraints w.r.t adv
            // from the input file
            opt_input::get_dobjcon_ds( mAbsDesVarVec, mOptCriteria, aGradObj_ds, aGradCon_ds );
        }

        // ---------------------------------------------------------------------

        void OptProb::assemble_grad(
                const Matrix< DDRMat >  & aGradObj_dCrit,
                const Matrix< DDRMat >  & aGradCon_dCrit,
                const Matrix< DDRMat >  & aGradObj_ds,
                const Matrix< DDRMat >  & aGradCon_ds )
        {
            // Assemble the implicit derivative of objective and constraints
            // w.r.t. adv i.e. the component that depends on the criteria which
            // in turn depends on the adv

            // Loop over number of advs
            for( uint i=0; i < mNumMasterAdvs; ++i )
            {
                // Loop over number of criteria
                for ( uint j=0; j < mNumCriteria; ++j )
                {
                    // dobj_ds    = dobj_ds       + dobj_dcrit          * dcrit_ds
                    mGradObj(i,0) = mGradObj(i,0) + aGradObj_dCrit(j,0) * mGradCrit(j,i);

                    // Loop over number of constraints
                    for ( uint k=0; k < mNumConstraints; ++k )
                    {
                        // dcon_ds    = dcon_ds       + dcon_dcrit          * dcon_ds
                        mGradCon(k,i) = mGradCon(k,i) + aGradCon_dCrit(k,j) * mGradCrit(j,i);
                    }
                }
            }

            // add the explicit components of the derivatives w.r.t. adv to the
            // gradients assembled above
            mGradObj = mGradObj + aGradObj_ds;
            mGradCon = mGradCon + aGradCon_ds;
        }

        // ---------------------------------------------------------------------

        void OptProb::scale_opt_sol()
        {
            // Need to decide on a framework for the scaling of the solution
        }

        // ---------------------------------------------------------------------

        void OptProb::update_opt_prob()
        {
            // Need to decide on a framework for the update of the problem
        }
    }
}
