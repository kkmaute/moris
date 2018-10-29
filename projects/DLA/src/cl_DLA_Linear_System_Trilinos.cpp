/*
 * cl_DLA_Linear_System_Trilinos.cpp
 *
 *  Created on: Dec 6, 2017
 *      Author: schmidt
 */
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"

#include "cl_DLA_Linear_System_Trilinos.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Enums.hpp"

using namespace moris;
using namespace dla;

Linear_System_Trilinos::Linear_System_Trilinos( Solver_Interface * aInput ) : moris::dla::Linear_Problem( aInput )
{
    if ( aInput->get_matrix_market_path() == NULL )
    {
        // Get number local dofs
        moris::uint aNumMyDofs = aInput->get_num_my_dofs();

        Matrix_Vector_Factory    tMatFactory;

        // create map object
        mMap = tMatFactory.create_map( aNumMyDofs,
                                       aInput->get_my_local_global_map(),
                                       aInput->get_constr_dof(),
                                       aInput->get_my_local_global_overlapping_map());      //FIXME

        // Build matrix
        mMat = tMatFactory.create_matrix( aInput, mMap );

        // Build RHS/LHS vector
        mVectorRHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );
        mFreeVectorLHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );

        mFullVectorLHS = tMatFactory.create_vector( aInput, mMap, VectorType::FULL_OVERLAPPING );

        mInput->build_graph( mMat );

        this->build_linear_system();
    }

    else
    {
        // Initialize communicator object
        Communicator_Epetra      tEpetraComm;

        // Build strings fot path to matrix, RHS and map
        char tMapString[500];      strcpy ( tMapString   , aInput->get_matrix_market_path() );   strcat( tMapString   , "map.mtx" );
        char tMatrixString[500];   strcpy ( tMatrixString, aInput->get_matrix_market_path() );   strcat( tMatrixString, "matrix.mtx" );
        char tVectorString[500];   strcpy ( tVectorString, aInput->get_matrix_market_path() );   strcat( tVectorString, "vector.mtx" );

        // Create map from matrix market file
        Epetra_Map * tMap;
        EpetraExt::MatrixMarketFileToMap( tMapString, *tEpetraComm.get_epetra_comm(), tMap );

        // Create matrix from matrix market file
        Epetra_CrsMatrix*   mMatFromMatrixMatket;
        //EpetraExt::MatrixMarketFileToCrsMatrix("/home/schmidt/matrix1.mtx", *tEpetraComm.get_epetra_comm(), mMatFromMatrixMatket);
        EpetraExt::MatrixMarketFileToCrsMatrix( tMatrixString, *tMap, mMatFromMatrixMatket );

        // Create RHS from matrix market file
        Epetra_MultiVector* mVecRHSFromMatrixMatket;
        EpetraExt::MatrixMarketFileToMultiVector( tVectorString , *tMap, mVecRHSFromMatrixMatket );

        Epetra_MultiVector* mVecLHSFromMatrixMatket;
        EpetraExt::MatrixMarketFileToMultiVector( tVectorString, *tMap, mVecLHSFromMatrixMatket );
        mVecLHSFromMatrixMatket->PutScalar(0.0);

        // Set matrix. solution vector and RHS
        mEpetraProblem.SetOperator( mMatFromMatrixMatket );
        mEpetraProblem.SetRHS( mVecRHSFromMatrixMatket );
        mEpetraProblem.SetLHS( mVecLHSFromMatrixMatket );
    }
}

//----------------------------------------------------------------------------------------

Linear_System_Trilinos::~Linear_System_Trilinos()
{
    delete( mMat );
    delete( mVectorRHS );
    delete( mFreeVectorLHS );
    delete( mMap );
}

//----------------------------------------------------------------------------------------
void Linear_System_Trilinos::assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector )
{
    mVectorRHS->vec_put_scalar( 0.0 );
    mMat->mat_put_scalar( 0.0 );

    mInput->fill_matrix_and_RHS( mMat, mVectorRHS, aFullSolutionVector);

//    mMat->print_matrix_to_screen();
//    std::cout<<*mVectorRHS->get_vector()<<std::endl;
}

//----------------------------------------------------------------------------------------
void Linear_System_Trilinos::assemble_residual( Dist_Vector * aFullSolutionVector )
{
    mVectorRHS->vec_put_scalar( 0.0 );

    mInput->assemble_RHS( mVectorRHS, aFullSolutionVector);

//    std::cout<<*mVectorRHS->get_vector()<<std::endl;
}

//----------------------------------------------------------------------------------------
void Linear_System_Trilinos::assemble_jacobian( Dist_Vector * aFullSolutionVector )
{
    mMat->mat_put_scalar( 0.0 );

    mInput->assemble_jacobian( mMat, aFullSolutionVector);

    //mMat->print_matrix_to_screen();
}

//----------------------------------------------------------------------------------------
void Linear_System_Trilinos::assemble_residual_and_jacobian( )
{
    mVectorRHS->vec_put_scalar( 0.0 );
    mMat->mat_put_scalar( 0.0 );

    mInput->fill_matrix_and_RHS( mMat, mVectorRHS);
}

//----------------------------------------------------------------------------------------
void Linear_System_Trilinos::build_linear_system()
 {
     // Set matrix. solution vector and RHS
     mEpetraProblem.SetOperator( mMat->get_matrix() );
     mEpetraProblem.SetRHS( mVectorRHS->get_vector() );
     mEpetraProblem.SetLHS( mFreeVectorLHS->get_vector() );

//     mMat->print_matrix_to_screen();
//     std::cout<<*mVectorRHS->get_vector()<<std::endl;
 }

//------------------------------------------------------------------------------------------
moris::sint Linear_System_Trilinos::solve_linear_system()
{
    moris::sint error = 0;
    // Get the linear system for the solver
    AztecOO Solver( mEpetraProblem );

    //Set solver options
    Solver.SetAztecOption( AZ_solver, AZ_gmres);
    Solver.SetAztecOption( AZ_precond, AZ_dom_decomp);
    Solver.SetAztecOption( AZ_diagnostics, 0);
    Solver.SetAztecOption( AZ_output, 0);

    //Solve
    error = Solver.Iterate( 200, 1E-8 );

    return error;
}

//------------------------------------------------------------------------------------------
void Linear_System_Trilinos::get_solution( Matrix< DDRMat > & LHSValues )
{
    mFreeVectorLHS->extract_copy( LHSValues );
}


