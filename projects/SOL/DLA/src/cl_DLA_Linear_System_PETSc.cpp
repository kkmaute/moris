/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_System_PETSc.cpp
 *
 */

#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_RowMatrix.h"

#include "cl_Vector_PETSc.hpp"
#include "cl_Vector_PETSc_Multi.hpp"
#include "cl_DLA_Linear_System_PETSc.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_SOL_Enums.hpp"
#include "cl_SOL_Warehouse.hpp"

#include <petsc.h>
#include <petscis.h>
#include <petscao.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <petsclog.h>

#include "cl_Stopwatch.hpp"    //CHR/src

using namespace moris;
using namespace dla;

Linear_System_PETSc::Linear_System_PETSc(
        Solver_Interface* aInput,
        const bool        aNotCreatedByNonLinSolver )
        : moris::dla::Linear_Problem( aInput )
        , mNotCreatedByNonLinearSolver( aNotCreatedByNonLinSolver )
{
    mTplType = sol::MapType::Petsc;
    if ( mNotCreatedByNonLinearSolver )
    {
        // Initialize petsc solvers
        PetscInitializeNoArguments();
    }

    if ( aInput->get_matrix_market_path() == NULL )
    {
        sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );

        // create map object
        mMap = tMatFactory.create_map(
                aInput->get_my_local_global_map(),
                aInput->get_constrained_Ids() );    // FIXME: should be full map?

        mMapFree = tMatFactory.create_map(
                aInput->get_my_local_global_map(),
                aInput->get_constrained_Ids() );    // FIXME

        // Build matrix
        mMat = tMatFactory.create_matrix( aInput, mMapFree );

        // Build RHS/LHS vector
        mFreeVectorLHS = tMatFactory.create_vector( aInput, mMapFree, 1 );


        mPointVectorRHS = tMatFactory.create_vector( aInput, mMapFree, 1 );
        mPointVectorLHS = tMatFactory.create_vector( aInput, mMapFree, 1 );

        mFullVectorLHS = tMatFactory.create_vector( aInput, mMap, 1 );    // FIXME: should be full map?

        // FIXME: graph not useful for petsc; needs to be done differently
        mSolverInterface->build_graph( mMat );
    }

    else
    {
        // MORIS_ERROR{ false, "Linear_System_PETSc::Linear_System_PETSc: Matrix market options not implemented in PETSc"};
    }
}

//------------------------------------------------------------------------------

Linear_System_PETSc::Linear_System_PETSc(
        Solver_Interface*   aInput,
        sol::SOL_Warehouse* aSolverWarehouse,
        sol::Dist_Map*      aFreeMap,
        sol::Dist_Map*      aFullMap,
        const bool          aNotCreatedByNonLinSolver )
        : moris::dla::Linear_Problem( aInput )
        , mNotCreatedByNonLinearSolver( aNotCreatedByNonLinSolver )
{
    mTplType         = sol::MapType::Petsc;
    mSolverWarehouse = aSolverWarehouse;

    if ( mNotCreatedByNonLinearSolver )
    {
        // Initialize petsc solvers
        PetscInitializeNoArguments();
    }

    sol::Matrix_Vector_Factory tMatFactory( sol::MapType::Petsc );

    // Build matrix
    mMat = tMatFactory.create_matrix( aInput, aFreeMap );

    uint tNumRHS = aInput->get_num_rhs();

    // Build RHS/LHS vector
    mFreeVectorLHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS );

    mPointVectorRHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS );
    mPointVectorLHS = tMatFactory.create_vector( aInput, aFreeMap, tNumRHS );

    mFullVectorLHS = tMatFactory.create_vector( aInput, aFullMap , tNumRHS);

    // start timer
    tic tTimer;

    mSolverInterface->build_graph( mMat, true );

    real tElapsedTime = tTimer.toc< moris::chronos::milliseconds >().wall;
    MORIS_LOG_INFO( "Building matrix graph on processor %u took %5.3f seconds.", (uint)par_rank(), (double)tElapsedTime / 1000 );
}

//------------------------------------------------------------------------------

Linear_System_PETSc::~Linear_System_PETSc()
{
    delete mMat;
    mMat = nullptr;

    delete mFreeVectorLHS;
    mFreeVectorLHS = nullptr;

    delete mFullVectorLHS;
    mFullVectorLHS = nullptr;

    delete mPointVectorLHS;
    mPointVectorLHS = nullptr;

    delete mPointVectorRHS;
    mPointVectorRHS = nullptr;

    //    mSolverInterface->delete_multigrid();

    // KSPDestroy( &mPetscProblem );
    //( &mpc );

    if ( mNotCreatedByNonLinearSolver == true )
    {
        // These calls are needed in order to delete the underlying Petsc maps before PetscFinalize
        delete mMap;
        delete mMapFree;

        PetscFinalize();
    }
}

//------------------------------------------------------------------------------

moris::sint
Linear_System_PETSc::solve_linear_system()
{
    // PetscInitializeNoArguments();
    KSP tPetscKSPProblem;
    PC  mpc;

    KSPCreate( PETSC_COMM_WORLD, &tPetscKSPProblem );
    KSPSetOperators( tPetscKSPProblem, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );
    KSPGetPC( tPetscKSPProblem, &mpc );

    // Build Preconditioner
    KSPGetPC( tPetscKSPProblem, &mpc );

    PCSetType( mpc, PCNONE );
    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
    PCFactorSetLevels( mpc, 0 );

    PetscInt maxits = 1000;
    KSPSetTolerances( tPetscKSPProblem, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
    KSPSetType( tPetscKSPProblem, KSPFGMRES );
    // KSPSetType(mPetscProblem,KSPPREONLY);
    KSPGMRESSetOrthogonalization( tPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPGMRESSetHapTol( tPetscKSPProblem, 1e-10 );
    KSPGMRESSetRestart( tPetscKSPProblem, 500 );

    KSPSetFromOptions( tPetscKSPProblem );

    // for debugging: print linear problem
    // KSPView( tPetscKSPProblem, PETSC_VIEWER_STDOUT_WORLD );

    // for debugging: print matrix and vectors
    // MatView( mMat->get_petsc_matrix(), PETSC_VIEWER_STDOUT_WORLD );
    // VecView( static_cast< Vector_PETSc* >( mPointVectorRHS )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );
    // VecView( static_cast< Vector_PETSc* >( mPointVectorLHS )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );

    Mat tRHSVecs = static_cast< MultiVector_PETSc* >( mPointVectorRHS )->get_petsc_vector();
    Mat tLHSVecs = static_cast< MultiVector_PETSc* >( mPointVectorLHS )->get_petsc_vector();

    for ( uint iNumRHS = 0; iNumRHS < mSolverInterface->get_num_rhs(); iNumRHS++ )
    {
        Vec tRHSVec, tLHSVec;
        
        MatDenseGetColumnVec( tRHSVecs, iNumRHS, &tRHSVec );
        MatDenseGetColumnVec( tLHSVecs, iNumRHS, &tLHSVec );

        // Create a new vector that has the same layout as the source vector
        VecAssemblyBegin( tRHSVec );
        VecAssemblyEnd( tRHSVec );
        VecAssemblyBegin( tLHSVec );
        VecAssemblyEnd( tLHSVec );


        KSPSolve( tPetscKSPProblem,tRHSVec,tLHSVec );


        MatDenseRestoreColumnVec( tRHSVecs, iNumRHS, &tRHSVec );
        MatDenseRestoreColumnVec( tLHSVecs, iNumRHS, &tLHSVec );
    }

    // VecView( static_cast< Vector_PETSc* >( mPointVectorLHS )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );
    KSPDestroy( &tPetscKSPProblem );

    return 0;
}

//------------------------------------------------------------------------------

void
Linear_System_PETSc::get_solution( Matrix< DDRMat >& LHSValues )
{
    mPointVectorLHS->extract_copy( LHSValues );
}

//------------------------------------------------------------------------------------------

void Linear_System_PETSc::construct_rhs_matrix()
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

    mSolverInterface->build_graph( mMassMat, true );
}