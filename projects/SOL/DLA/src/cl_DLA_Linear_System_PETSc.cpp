/*
 * cl_DLA_Linear_System_PETSc.cpp
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

#include "cl_DLA_Linear_System_PETSc.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Enums.hpp"

#include <petsc.h>
#include <petscis.h>
#include <petscao.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <petsclog.h>

using namespace moris;
using namespace dla;

Linear_System_PETSc::Linear_System_PETSc(       Solver_Interface * aInput,
                                          const bool               aNotCreatedByNonLinSolver) : moris::dla::Linear_Problem( aInput )
{
    mNotCreatedByNonLinearSolver = aNotCreatedByNonLinSolver;

    // Initialize petsc solvers
    PetscInitializeNoArguments();

    if ( aInput->get_matrix_market_path() == NULL )
    {
        // Initialize petsc solvers
        //PetscInitializeNoArguments();

        // Get number local dofs
        moris::uint aNumMyDofs = aInput->get_num_my_dofs();

        Matrix_Vector_Factory tMatFactory( MapType::Petsc );

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
        //MORIS_ERROR{ false, "Linear_System_PETSc::Linear_System_PETSc: Matrix market options not implemented in PETSc"};
    }
}

//----------------------------------------------------------------------------------------

Linear_System_PETSc::~Linear_System_PETSc()
{
    delete( mMat );
    delete( mVectorRHS );
    delete( mFreeVectorLHS );
    delete( mFullVectorLHS );
    delete( mMap );

    //KSPDestroy( &mPetscProblem );
    //( &mpc );

    if ( mNotCreatedByNonLinearSolver == true )
    {
        PetscFinalize();
    }
}

//----------------------------------------------------------------------------------------
void Linear_System_PETSc::build_linear_system()
 {
     // Set matrix. solution vector and RHS
     //KSPSetOperators( tPetscKSPProblem, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );

//     mMat->print_matrix_to_screen();
//     std::cout<<*mVectorRHS->get_vector()<<std::endl;
 }

//------------------------------------------------------------------------------------------
moris::sint Linear_System_PETSc::solve_linear_system()
{
    //PetscInitializeNoArguments();
    KSP tPetscKSPProblem;
    PC mpc;

    KSPCreate( PETSC_COMM_WORLD, &tPetscKSPProblem );
    KSPSetOperators( tPetscKSPProblem, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );
    //KSPView ( tPetscKSPProblem, PETSC_VIEWER_STDOUT_WORLD);
    KSPGetPC( tPetscKSPProblem, &mpc );

    // Build Preconditioner
    KSPGetPC( tPetscKSPProblem, &mpc );

    PCSetType( mpc, PCNONE );
    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
    PCFactorSetLevels( mpc, 0 );

    PetscInt maxits=1000;
    KSPSetTolerances( tPetscKSPProblem, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
    KSPSetType( tPetscKSPProblem, KSPFGMRES );
    //KSPSetType(mPetscProblem,KSPPREONLY);
    KSPGMRESSetOrthogonalization( tPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPGMRESSetHapTol( tPetscKSPProblem, 1e-10 );
    KSPGMRESSetRestart( tPetscKSPProblem, 500 );

    KSPSetFromOptions( tPetscKSPProblem );

    KSPSolve( tPetscKSPProblem, mVectorRHS->get_petsc_vector(), mFreeVectorLHS->get_petsc_vector() );

    KSPDestroy( &tPetscKSPProblem );

    return 0;
}

//------------------------------------------------------------------------------------------
void Linear_System_PETSc::get_solution( Matrix< DDRMat > & LHSValues )
{
    //VecGetArray (tSolution, &  LHSValues.data());

    moris::sint tVecLocSize;
    VecGetLocalSize( mFreeVectorLHS->get_petsc_vector(), &tVecLocSize );

    // FIXME replace with VecGetArray()
    moris::Matrix< DDSMat > tVal ( tVecLocSize, 1 , 0 );
    LHSValues.set_size( tVecLocSize, 1 );

    //----------------------------------------------------------------------------------------
    // Get list containing the number of owned adofs of each processor
    Matrix< DDUMat > tNumOwnedList;
    comm_gather_and_broadcast( tVecLocSize, tNumOwnedList );

    Matrix< DDUMat > tOwnedOffsetList( tNumOwnedList.length(), 1, 0 );

    // Loop over all entries to create the offsets. Starting with 1
    for ( moris::uint Ij = 1; Ij < tOwnedOffsetList.length(); Ij++ )
    {
        // Add the number of owned adofs of the previous processor to the offset of the previous processor
        tOwnedOffsetList( Ij, 0 ) = tOwnedOffsetList( Ij-1, 0 ) + tNumOwnedList( Ij-1, 0 );
    }
    //-------------------------------------------------------------------------------------
    for ( moris::sint Ik=0; Ik< tVecLocSize; Ik++ )
    {
        tVal( Ik, 0 ) = tOwnedOffsetList( par_rank(), 0)+Ik;
    }

    VecGetValues( mFreeVectorLHS->get_petsc_vector(), tVecLocSize, tVal.data(), LHSValues.data() );
}


