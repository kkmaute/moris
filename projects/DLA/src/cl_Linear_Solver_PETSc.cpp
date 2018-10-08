/*
 * cl_Linear_Solver_PETSc.cpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#include "cl_Linear_Solver_PETSc.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

using namespace moris;

moris::Linear_Solver_PETSc::Linear_Solver_PETSc( moris::Solver_Interface * aInput ) : moris::Linear_Solver( aInput )
{
    // Initialize petsc solvers
    PetscInitializeNoArguments();
    KSPCreate( PETSC_COMM_WORLD, &mksp );
    KSPGetPC( mksp, &mpc );

    moris::uint aNumMyDofs = aInput->get_num_my_dofs();

    moris::Matrix_Vector_Factory      tMatFactory;

    // create map object
    mMap = tMatFactory.create_map( aNumMyDofs,
                                   aInput->get_my_local_global_map(),
                                   aInput->get_constr_dof(),
                                   aInput->get_my_local_global_map());

    // Build matrix
    mMat = tMatFactory.create_matrix( aInput, mMap );

    // Build RHS/LHS vector
    mVectorRHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );
    mVectorLHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );

    moris::Model_Solver_Interface tLinProblem( this, aInput, mMat, mVectorRHS );
}

//moris::Linear_Solver_PETSc::Linear_Solver_PETSc( Mat          aPETScMat,
//                                          Vec          aPETScVector_x,
//                                          Vec          aPETScVector_b)
//{
//
//    // Build linear system
//    KSPCreate( PETSC_COMM_WORLD, &mksp );
//    KSPSetOperators( mksp, aPETScMat, aPETScMat );
//
//    // Build Preconditioner
//    KSPGetPC( mksp, &mpc );
//
//    PCSetType( mpc, PCNONE );
//    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
//    PCFactorSetLevels( mpc, 0 );
//
//    PetscInt maxits=1000;
//    KSPSetTolerances( mksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
//    KSPSetType( mksp, KSPFGMRES );
//    //KSPSetType(mksp,KSPPREONLY);
//    KSPGMRESSetOrthogonalization( mksp, KSPGMRESModifiedGramSchmidtOrthogonalization );
//    KSPGMRESSetHapTol( mksp, 1e-10 );
//    KSPGMRESSetRestart( mksp, 500 );
//
//    KSPSetFromOptions( mksp );
//
//    KSPSolve( mksp, aPETScVector_b, aPETScVector_x );
//}

void moris::Linear_Solver_PETSc::build_linear_system()
{
    // build linear system
    KSPSetOperators( mksp, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );
}

moris::sint moris::Linear_Solver_PETSc::solve_linear_system()
{
    // set Petsc preconditioner
    PCSetType( mpc, PCNONE );
    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
    PCFactorSetLevels( mpc, 0 );

    PetscInt maxits=1000;
    KSPSetTolerances( mksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
    KSPSetType(mksp,KSPFGMRES);
    //KSPSetType(mksp,KSPPREONLY);
    KSPGMRESSetOrthogonalization( mksp, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPGMRESSetHapTol( mksp, 1e-10 );
    KSPGMRESSetRestart( mksp, 500 );

    KSPSetFromOptions( mksp );

    KSPSolve( mksp, mVectorRHS->get_petsc_vector(), mVectorLHS->get_petsc_vector() );

    return 0;
}

moris::Linear_Solver_PETSc::~Linear_Solver_PETSc()
{
    KSPDestroy( &mksp );
    //PCDestroy( &mpc );

    delete( mMat );
    delete( mVectorRHS );
    delete( mVectorLHS );
    delete( mMap );

    PetscFinalize();
}

void moris::Linear_Solver_PETSc::get_solution( moris::Matrix< DDRMat > & LHSValues )
{
    Vec tSolution;

    KSPGetSolution( mksp, &tSolution );

    //VecGetArray (tSolution, &  LHSValues.data());

    // FIXME replace with VecGetArray()
    moris::Matrix< DDSMat > tVal ( LHSValues.length(), 1 );

    for ( moris::uint Ik=0; Ik< LHSValues.length(); Ik++ )
    {
        tVal( Ik, 0 ) = Ik;
    }

    //VecView ( tSolution, PETSC_VIEWER_STDOUT_WORLD);

    VecGetValues ( tSolution, LHSValues.length(), tVal.data() , LHSValues.data() );

    //VecDestroy( &tSolution );

}
