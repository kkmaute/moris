/*
 * cl_DLA_Linear_Problem.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_Vector.hpp"
#include "cl_Sparse_Matrix.hpp"

namespace moris
{
namespace dla
{
    Dist_Vector * Linear_Problem::get_full_solver_LHS()
    {
        mFullVectorLHS->vec_put_scalar( 0.0 );

        mFullVectorLHS->import_local_to_global( *mFreeVectorLHS );

        return mFullVectorLHS;
    }

//----------------------------------------------------------------------------------------
    void Linear_Problem::set_free_solver_LHS( Dist_Vector * aFullSolVector)
    {
        mFreeVectorLHS->import_local_to_global( *aFullSolVector );
    }

//----------------------------------------------------------------------------------------
    void Linear_Problem::assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector )
    {
        mVectorRHS->vec_put_scalar( 0.0 );
        mMat->mat_put_scalar( 0.0 );

        mInput->fill_matrix_and_RHS( mMat, mVectorRHS, aFullSolutionVector);

        //mMat->print_matrix_to_screen();
        //std::cout<<*mVectorRHS->get_vector()<<std::endl;
    }

//----------------------------------------------------------------------------------------
    void Linear_Problem::assemble_residual( Dist_Vector * aFullSolutionVector )
    {
        mVectorRHS->vec_put_scalar( 0.0 );

        mInput->assemble_RHS( mVectorRHS, aFullSolutionVector);

        //std::cout<<*mVectorRHS->get_vector()<<std::endl;
    }

//----------------------------------------------------------------------------------------
    void Linear_Problem::assemble_jacobian( Dist_Vector * aFullSolutionVector )
    {
        mMat->mat_put_scalar( 0.0 );

        mInput->assemble_jacobian( mMat, aFullSolutionVector);

        //mMat->print_matrix_to_screen();
    }

//----------------------------------------------------------------------------------------
    void Linear_Problem::assemble_residual_and_jacobian( )
    {
        mVectorRHS->vec_put_scalar( 0.0 );
        mMat->mat_put_scalar( 0.0 );

        mInput->fill_matrix_and_RHS( mMat, mVectorRHS);
    }

}
}

