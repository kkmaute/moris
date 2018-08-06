/*
 * cl_Linear_Solver_PETSc.cpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#include "cl_Linear_Solver_PETSc.hpp"
#include "cl_Matrix_Vector_Factory.hpp"
//#include "cl_Solver_Input.hpp"
/*
moris::Linear_Solver_PETSc::Linear_Solver_PETSc( moris::Solver_Input * aInput ) : moris::Linear_Solver()
{
    PetscInitializeNoArguments();
    KSPCreate( PETSC_COMM_WORLD, &mksp );
    KSPGetPC( mksp, &mpc );

    moris::uint aNumMyDofs = aInput->get_num_my_dofs();

    moris::Matrix_Vector_Factory      tMatFactory;

    // create map object
    mMap = tMatFactory.create_map( aNumMyDofs,
                                   aInput->get_my_local_global_map(),
                                   aInput->get_constr_dof() );

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
    KSPSetOperators( mksp, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );
}

void moris::Linear_Solver_PETSc::solve_linear_system()
{
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

void moris::Linear_Solver_PETSc::get_solution( moris::Mat< moris::real > & LHSValues )
{
    Vec tSolution;

    KSPGetSolution( mksp, &tSolution );

    //VecGetArray (tSolution, & mem_pointer( LHSValues ));

    // FIXME replace with VecGetArray()
    moris::Mat < int > tVal ( LHSValues.length(), 1 );

    for ( moris::uint Ik=0; Ik< LHSValues.length(); Ik++ )
    {
        tVal( Ik, 0 ) = Ik;
    }

    //VecView ( tSolution, PETSC_VIEWER_STDOUT_WORLD);

    VecGetValues ( tSolution, LHSValues.length(), mem_pointer( tVal ) ,mem_pointer( LHSValues ) );

    //VecDestroy( &tSolution );

}*/
