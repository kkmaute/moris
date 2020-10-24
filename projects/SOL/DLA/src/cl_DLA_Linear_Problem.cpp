/*
 * cl_DLA_Linear_Problem.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_Stopwatch.hpp" //CHR/src

// detailed logging package
#include "cl_Tracer.hpp"
#include "cl_Tracer_Enums.hpp"

namespace moris
{
    namespace dla
    {
        sol::Dist_Vector * Linear_Problem::get_full_solver_LHS()
        {
            // zero out full LHS vec
            mFullVectorLHS->vec_put_scalar( 0.0 );

            // Import free LHS to full LHS
            mFullVectorLHS->import_local_to_global( *mFreeVectorLHS );

            return mFullVectorLHS;
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::set_free_solver_LHS( sol::Dist_Vector * aFullSolVector)
        {
            mFreeVectorLHS->import_local_to_global( *aFullSolVector );
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual_and_jacobian( sol::Dist_Vector * aFullSolutionVector )
        {
            // zero out RHS
            mVectorRHS->vec_put_scalar( 0.0 );

            // zero out matrix
            mMat->mat_put_scalar( 0.0 );

            //
            mSolverInterface->fill_matrix_and_RHS( mMat, mVectorRHS, aFullSolutionVector );
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual()
        {
            Tracer tTracer(EntityBase::LinearProblem, EntityType::NoType, EntityAction::AssembleResidual);

            // Zero out RHS
            mVectorRHS->vec_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            // assemble RHS
            mSolverInterface->assemble_RHS( mVectorRHS );

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mVectorRHS );

            if( !mSolverInterface->get_is_forward_analysis() )
            {
                this->compute_residual_for_adjoint_solve();
            }


            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
            MORIS_LOG_INFO( " Assembly of residual on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);

            //std::cout<<"first RHS"<<std::endl;
            //mVectorRHS->print();
        }

        void Linear_Problem::assemble_staggered_residual_contribution()
        {
            mSolverInterface->assemble_staggerd_RHS_contribution( mVectorRHS );
            //std::cout<<"second RHS"<<std::endl;
            //mVectorRHS->print();
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_jacobian()
        {
            Tracer tTracer(EntityBase::LinearProblem, EntityType::NoType, EntityAction::AssembleJacobian);

            mMat->mat_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            // assemble Jacobian
            mSolverInterface->assemble_jacobian( mMat);

            // stop timer
            real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;

            MORIS_LOG_INFO( " Assembly of Jacobian on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);
        }

        //----------------------------------------------------------------------------------------
        void Linear_Problem::assemble_residual_and_jacobian( )
        {
            mVectorRHS->vec_put_scalar( 0.0 );
            mMat->mat_put_scalar( 0.0 );

            // start timer
            tic tTimer;

            mSolverInterface->fill_matrix_and_RHS( mMat, mVectorRHS );

            mSolverInterface->assemble_additional_DqDs_RHS_contribution( mVectorRHS );

            if( !mSolverInterface->get_is_forward_analysis() )
            {
                this->compute_residual_for_adjoint_solve();
            }

            // stop timer
            //real tElapsedTime = tTimer.toc<moris::chronos::milliseconds>().wall;
            //MORIS_LOG_INFO( "Assembly of Residual and Jacobian on processor %u took %5.3f seconds.", ( uint ) par_rank(), ( double ) tElapsedTime / 1000);
        }

        //----------------------------------------------------------------------------------------

        void Linear_Problem::compute_residual_for_adjoint_solve()
        {
            sol::Matrix_Vector_Factory tMatFactory( mTplType );

            // create vector for jacobian times sol vec
            sol::Dist_Vector * tMatTimesSolVec = tMatFactory.create_vector(
                          mSolverInterface,
                          mVectorRHS->get_map(),
                          mSolverInterface->get_num_rhs() );

            // multiply jacobian with previous solution vector
            mMat->mat_vec_product( *mFreeVectorLHS, *tMatTimesSolVec, false );

            // add contribution to RHS
            mVectorRHS->vec_plus_vec( 1.0, *tMatTimesSolVec, 1.0 );

            delete tMatTimesSolVec;
        }
    }
}

