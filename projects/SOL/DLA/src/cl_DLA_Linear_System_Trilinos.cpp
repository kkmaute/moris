/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_System_Trilinos.cpp
 *
 */

#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"

#include "cl_DLA_Linear_System_Trilinos.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Enums.hpp"

#include "cl_Stopwatch.hpp"    //CHR/src

// detailed logging package
#include "cl_Tracer.hpp"

using namespace moris;
using namespace dla;

Linear_System_Trilinos::Linear_System_Trilinos( Solver_Interface* aInput )
        : moris::dla::Linear_Problem( aInput )
{
    mTplType = sol::MapType::Epetra;
    if ( aInput->get_matrix_market_path() == NULL )
    {
        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Epetra );

        // create map object
        mMapFree = tMatFactory.create_map(
                aInput->get_my_local_global_map(),
                aInput->get_constrained_Ids() );    // FIXME

        mMap = tMatFactory.create_full_map(
                aInput->get_my_local_global_map(),
                aInput->get_my_local_global_overlapping_map() );    // FIXME

        mMapFree->build_dof_translator( aInput->get_my_local_global_overlapping_map(), false );

        // Build matrix
        mMat = tMatFactory.create_matrix( aInput, mMapFree, true, true );

        uint tNumRHS = aInput->get_num_rhs();

        // Build RHS/LHS vector
        mFreeVectorLHS = tMatFactory.create_vector( aInput, mMapFree, tNumRHS );

        mPointVectorRHS = tMatFactory.create_vector( aInput, mMapFree, tNumRHS, true );
        mPointVectorLHS = tMatFactory.create_vector( aInput, mMapFree, tNumRHS, true );

        mFullVectorLHS = tMatFactory.create_vector( aInput, mMap, tNumRHS );

        // start timer
        tic tTimer;

        mSolverInterface->build_graph( mMat );

        real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;
        MORIS_LOG_INFO( "Building matrix graph on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );
    }

    else
    {
        //        // Initialize communicator object
        //        Communicator_Epetra      tEpetraComm;
        //
        //        // Build strings fot path to matrix, RHS and map
        //        char tMapString[500];      strcpy ( tMapString   , aInput->get_matrix_market_path() );   strcat( tMapString   , "map.mtx" );
        //        char tMatrixString[500];   strcpy ( tMatrixString, aInput->get_matrix_market_path() );   strcat( tMatrixString, "matrix.mtx" );
        //        char tVectorString[500];   strcpy ( tVectorString, aInput->get_matrix_market_path() );   strcat( tVectorString, "vector.mtx" );
        //
        //        // Create map from matrix market file
        //        Epetra_Map * tMap;
        //        EpetraExt::MatrixMarketFileToMap( tMapString, *tEpetraComm.get_epetra_comm(), tMap );
        //
        //        // Create matrix from matrix market file
        //        Epetra_CrsMatrix*   mMatFromMatrixMatket;
        //        //EpetraExt::MatrixMarketFileToCrsMatrix("/home/schmidt/matrix1.mtx", *tEpetraComm.get_epetra_comm(), mMatFromMatrixMatket);
        //        EpetraExt::MatrixMarketFileToCrsMatrix( tMatrixString, *tMap, mMatFromMatrixMatket );
        //
        //        // Create RHS from matrix market file
        //        Epetra_MultiVector* mVecRHSFromMatrixMatket;
        //        EpetraExt::MatrixMarketFileToMultiVector( tVectorString , *tMap, mVecRHSFromMatrixMatket );
        //
        //        Epetra_MultiVector* mVecLHSFromMatrixMatket;
        //        EpetraExt::MatrixMarketFileToMultiVector( tVectorString, *tMap, mVecLHSFromMatrixMatket );
        //        mVecLHSFromMatrixMatket->PutScalar(0.0);
        //
        //        mMat->get_matrix()             = mMatFromMatrixMatket;
        //		mVectorRHS->get_vector()       = mVecRHSFromMatrixMatket;
        //		mFreeVectorLHS->get_vector()   = mVecLHSFromMatrixMatket;
    }
}

//----------------------------------------------------------------------------------------
Linear_System_Trilinos::Linear_System_Trilinos(
        Solver_Interface*   aInput,
        sol::SOL_Warehouse* aSolverWarehouse,
        sol::Dist_Map*      aFreeMap,
        sol::Dist_Map*      aFullMap )
        : moris::dla::Linear_Problem( aInput )
{
    mTplType         = sol::MapType::Epetra;
    mSolverWarehouse = aSolverWarehouse;
    sol::Matrix_Vector_Factory tMatFactory( mTplType );

    aFreeMap->build_dof_translator( aInput->get_my_local_global_overlapping_map(), false );

    // Build matrix
    mMat = tMatFactory.create_matrix( aInput, aFreeMap, true, true );

    uint tNumRHS = aInput->get_num_rhs();

    // Build RHS/LHS vector
    mFreeVectorLHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS );

    mPointVectorRHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS, true );
    mPointVectorLHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS, true );

    mFullVectorLHS = tMatFactory.create_vector( aInput, aFullMap, tNumRHS );

    // start timer
    tic tTimer;

    mSolverInterface->build_graph( mMat );

    real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;
    MORIS_LOG_INFO( "Building matrix graph on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );
}

//----------------------------------------------------------------------------------------

Linear_System_Trilinos::~Linear_System_Trilinos()
{
    delete mMat;
    mMat = nullptr;

    delete mMassMat;
    mMassMat = nullptr;

    delete mFreeVectorLHS;
    mFreeVectorLHS = nullptr;

    delete mFullVectorLHS;
    mFullVectorLHS = nullptr;

    delete mPointVectorLHS;
    mPointVectorLHS = nullptr;

    delete mPointVectorRHS;
    mPointVectorRHS = nullptr;

    delete mMap;
    delete mMapFree;
}

//------------------------------------------------------------------------------------------
moris::sint
Linear_System_Trilinos::solve_linear_system()
{
    moris::sint error = 0;
    // Get the linear system for the solver

    Epetra_LinearProblem tEpetraProblem;
    tEpetraProblem.SetOperator( mMat->get_matrix() );
    tEpetraProblem.SetRHS( static_cast< Vector_Epetra* >( mPointVectorRHS )->get_epetra_vector() );
    tEpetraProblem.SetLHS( static_cast< Vector_Epetra* >( mPointVectorLHS )->get_epetra_vector() );

    AztecOO Solver( tEpetraProblem );

    // Set solver options
    Solver.SetAztecOption( AZ_solver, AZ_gmres );
    Solver.SetAztecOption( AZ_precond, AZ_dom_decomp );
    Solver.SetAztecOption( AZ_diagnostics, 0 );
    Solver.SetAztecOption( AZ_output, 0 );

    // Solve
    error = Solver.Iterate( 200, 1E-8 );

    return error;
}

//------------------------------------------------------------------------------------------
void Linear_System_Trilinos::get_solution( Matrix< DDRMat >& LHSValues )
{
    mPointVectorLHS->extract_copy( LHSValues );
}

//------------------------------------------------------------------------------------------
void Linear_System_Trilinos::construct_rhs_matrix()
{
       //delete the previous mass matrix if it exits
    if ( mMassMat != nullptr )
    {
        delete mMassMat;
        mMassMat = nullptr;
    }
    
    // use copy constructor to create mass matrix
    sol::Matrix_Vector_Factory tMatFactory( mSolverWarehouse->get_tpl_type() );

    // Build matrix
    mMassMat = tMatFactory.create_matrix( mSolverInterface, mMat->get_map(), true, true );

    mSolverInterface->build_graph( mMassMat );
}
