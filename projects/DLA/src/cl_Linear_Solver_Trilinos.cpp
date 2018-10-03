/*
 * LinearSolverTrilinos.cpp
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

//#include "cl_MSI_Model_Solver_Interface.hpp"

#include "cl_Linear_Solver_Trilinos.hpp"
#include "cl_Solver_Input.hpp"
#include "cl_DistLinAlg_Enums.hpp"

using namespace moris;

Linear_Solver_Trilinos::Linear_Solver_Trilinos( Solver_Input *  aInput ) : moris::Linear_Solver( aInput ),
                                                                           mEpetraProblem()
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
        mVectorLHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );

        mVectorLHSOverlapping = tMatFactory.create_vector( aInput, mMap, VectorType::FULL_OVERLAPPING );

        Model_Solver_Interface tLinProblem;

        tLinProblem.build_graph( this->get_solver_input(), mMat );

        this->build_linear_system();

        //Model_Solver_Interface tLinProblem( this, aInput, mMat, mVectorRHS );
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

Linear_Solver_Trilinos::~Linear_Solver_Trilinos()
{
    delete( mMat );
    delete( mVectorRHS );
    delete( mVectorLHS );
    delete( mVectorLHSOverlapping );
    delete( mMap );
}

//----------------------------------------------------------------------------------------
void Linear_Solver_Trilinos::assemble_residual_and_jacobian( Dist_Vector * aFullSolutionVector )
{
    mVectorRHS->vec_put_scalar( 0.0 );
    mMat->mat_put_scalar( 0.0 );
    //Model_Solver_Interface tLinProblem( this, this->get_solver_input(), mMat, mVectorRHS );
    Model_Solver_Interface tLinProblem;
    tLinProblem.fill_matrix_and_RHS( this, this->get_solver_input(), mMat, mVectorRHS, aFullSolutionVector);
}

//----------------------------------------------------------------------------------------
void Linear_Solver_Trilinos::assemble_residual_and_jacobian( )
{
    mVectorRHS->vec_put_scalar( 0.0 );
    mMat->mat_put_scalar( 0.0 );
    //Model_Solver_Interface tLinProblem( this, this->get_solver_input(), mMat, mVectorRHS );
    Model_Solver_Interface tLinProblem;
    tLinProblem.fill_matrix_and_RHS(this, this->get_solver_input(), mMat, mVectorRHS);
}

//----------------------------------------------------------------------------------------
void Linear_Solver_Trilinos::build_linear_system()
 {
     // Set matrix. solution vector and RHS
     mEpetraProblem.SetOperator( mMat->get_matrix() );
     mEpetraProblem.SetRHS( mVectorRHS->get_vector() );
     mEpetraProblem.SetLHS( mVectorLHS->get_vector() );

//     mMat->print_matrix_to_screen();
//     std::cout<<*mVectorRHS->get_vector()<<std::endl;
 }

//------------------------------------------------------------------------------------------
moris::sint Linear_Solver_Trilinos::solve_linear_system()
{
    moris::sint error = 0;
    // Get the linear system for the solver
    AztecOO Solver( mEpetraProblem );

    //Set solver options
    Solver.SetAztecOption( AZ_solver, AZ_gmres);
    Solver.SetAztecOption( AZ_precond, AZ_dom_decomp);

    //Solve
    error = Solver.Iterate( 200, 1E-8 );

    //std::cout << "Solver performed " << Solver.NumIters()  << "iterations.\n";
    //std::cout << "Norm of the true residual = " << Solver.TrueResidual() << std::endl;
    return error;
}

//------------------------------------------------------------------------------------------
void Linear_Solver_Trilinos::get_solution( Matrix< DDRMat > & LHSValues )
{
    mVectorLHS->extract_copy( LHSValues );
}

//-------------------------------------------------------------------------------------------
void Linear_Solver_Trilinos::solve_eigenvalues()
{
}

