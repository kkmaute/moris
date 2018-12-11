/*
 * cl_DLA_Linear_Solver_PETSc.cpp
 *
 *  Created on: Mar 25, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Solver_PETSc.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_PETSc::Linear_Solver_PETSc()
{
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::Linear_Solver_PETSc( std::shared_ptr< Linear_Problem > aLinearSystem )
{
    mLinearSystem = aLinearSystem;

    //FIXME add rest
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::~Linear_Solver_PETSc()
{

}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_linear_problem( std::shared_ptr< Linear_Problem > aLinearSystem )
{
    mLinearSystem = aLinearSystem;
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_parameters()
{

}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system( )
{
    this->set_solver_internal_parameters();

    // Build linear system
    KSPCreate( PETSC_COMM_WORLD, &mPetscKSPProblem );
    KSPSetOperators( mPetscKSPProblem, mLinearSystem->get_matrix()->get_petsc_matrix(), mLinearSystem->get_matrix()->get_petsc_matrix() );

    PC mpc;

    // Build Preconditioner
    KSPGetPC( mPetscKSPProblem, &mpc );

    PCSetType( mpc, PCNONE );
    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
    PCFactorSetLevels( mpc, 0 );

    PetscInt maxits=1000;
    KSPSetTolerances( mPetscKSPProblem, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
    KSPSetType( mPetscKSPProblem, KSPFGMRES );
    //KSPSetType(mPetscKSPProblem,KSPPREONLY);
    KSPGMRESSetOrthogonalization( mPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPGMRESSetHapTol( mPetscKSPProblem, 1e-10 );
    KSPGMRESSetRestart( mPetscKSPProblem, 500 );

    KSPSetFromOptions( mPetscKSPProblem );

    KSPSolve( mPetscKSPProblem, mLinearSystem->get_solver_RHS()->get_petsc_vector(), mLinearSystem->get_free_solver_LHS()->get_petsc_vector() );


    return 0;
}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system(       std::shared_ptr< Linear_Problem > aLinearSystem,
                                                      const moris::sint                       aIter )
{
    mLinearSystem = aLinearSystem;

    this->set_solver_internal_parameters();

    // Build linear system
    PC mpc;

    KSPCreate( PETSC_COMM_WORLD, &mPetscKSPProblem );

    KSPSetOperators( mPetscKSPProblem, mLinearSystem->get_matrix()->get_petsc_matrix(), mLinearSystem->get_matrix()->get_petsc_matrix() );


    // Build Preconditioner
    KSPGetPC( mPetscKSPProblem, &mpc );

    PCSetType( mpc, PCNONE );
    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
    PCFactorSetLevels( mpc, 0 );

    PetscInt maxits=1000;
    KSPSetTolerances( mPetscKSPProblem, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
    KSPSetType( mPetscKSPProblem, KSPFGMRES );
    //KSPSetType(mksp,KSPPREONLY);
    KSPGMRESSetOrthogonalization( mPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPGMRESSetHapTol( mPetscKSPProblem, 1e-10 );
    KSPGMRESSetRestart( mPetscKSPProblem, 500 );

    KSPSetFromOptions( mPetscKSPProblem );

    KSPSolve( mPetscKSPProblem, mLinearSystem->get_solver_RHS()->get_petsc_vector(), mLinearSystem->get_free_solver_LHS()->get_petsc_vector() );

    return 0;
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_internal_parameters()
{

}

//moris::Linear_Solver_PETSc::Linear_Solver_PETSc( moris::Solver_Interface * aInput ) //: moris::Linear_Solver( aInput )
//{
//    // Initialize petsc solvers
//    PetscInitializeNoArguments();
//    KSPCreate( PETSC_COMM_WORLD, &mksp );
//    KSPGetPC( mksp, &mpc );
//
//    moris::uint aNumMyDofs = aInput->get_num_my_dofs();
//
//    moris::Matrix_Vector_Factory      tMatFactory;
//
//    // create map object
//    mMap = tMatFactory.create_map( aNumMyDofs,
//                                   aInput->get_my_local_global_map(),
//                                   aInput->get_constr_dof(),
//                                   aInput->get_my_local_global_map());
//
//    // Build matrix
//    mMat = tMatFactory.create_matrix( aInput, mMap );
//
//    // Build RHS/LHS vector
//    mVectorRHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );
//    mFreeVectorLHS = tMatFactory.create_vector( aInput, mMap, VectorType::FREE );
//
//    mInput->build_graph( mMat );
//}
//
////moris::Linear_Solver_PETSc::Linear_Solver_PETSc( Mat          aPETScMat,
////                                          Vec          aPETScVector_x,
////                                          Vec          aPETScVector_b)
////{
////
////    // Build linear system
////    KSPCreate( PETSC_COMM_WORLD, &mksp );
////    KSPSetOperators( mksp, aPETScMat, aPETScMat );
////
////    // Build Preconditioner
////    KSPGetPC( mksp, &mpc );
////
////    PCSetType( mpc, PCNONE );
////    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
////    PCFactorSetLevels( mpc, 0 );
////
////    PetscInt maxits=1000;
////    KSPSetTolerances( mksp, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
////    KSPSetType( mksp, KSPFGMRES );
////    //KSPSetType(mksp,KSPPREONLY);
////    KSPGMRESSetOrthogonalization( mksp, KSPGMRESModifiedGramSchmidtOrthogonalization );
////    KSPGMRESSetHapTol( mksp, 1e-10 );
////    KSPGMRESSetRestart( mksp, 500 );
////
////    KSPSetFromOptions( mksp );
////
////    KSPSolve( mksp, aPETScVector_b, aPETScVector_x );
////}
//
//void moris::Linear_Solver_PETSc::build_linear_system()
//{
//    // build linear system
//    KSPSetOperators( mksp, mMat->get_petsc_matrix(), mMat->get_petsc_matrix() );
//}
//
//moris::sint moris::Linear_Solver_PETSc::solve_linear_system()
//{
//    // set Petsc preconditioner
//    PCSetType( mpc, PCNONE );
//    PCFactorSetDropTolerance( mpc, 1e-6, PETSC_DEFAULT, PETSC_DEFAULT );
//    PCFactorSetLevels( mpc, 0 );
//
//    PetscInt maxits=1000;
//    KSPSetTolerances( mksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, maxits );
//    KSPSetType(mksp,KSPFGMRES);
//    //KSPSetType(mksp,KSPPREONLY);
//    KSPGMRESSetOrthogonalization( mksp, KSPGMRESModifiedGramSchmidtOrthogonalization );
//    KSPGMRESSetHapTol( mksp, 1e-10 );
//    KSPGMRESSetRestart( mksp, 500 );
//
//    KSPSetFromOptions( mksp );
//
//    KSPSolve( mksp, mVectorRHS->get_petsc_vector(), mFreeVectorLHS->get_petsc_vector() );
//
//    return 0;
//}
//
//moris::Linear_Solver_PETSc::~Linear_Solver_PETSc()
//{
//    KSPDestroy( &mksp );
//    //PCDestroy( &mpc );
//
//    delete( mMat );
//    delete( mVectorRHS );
//    delete( mFreeVectorLHS );
//    delete( mMap );
//
//    PetscFinalize();
//}
//
//void moris::Linear_Solver_PETSc::get_solution( moris::Matrix< DDRMat > & LHSValues )
//{
//    Vec tSolution;
//
//    KSPGetSolution( mksp, &tSolution );
//
//    //VecGetArray (tSolution, &  LHSValues.data());
//
//    // FIXME replace with VecGetArray()
//    moris::Matrix< DDSMat > tVal ( LHSValues.length(), 1 );
//
//    for ( moris::uint Ik=0; Ik< LHSValues.length(); Ik++ )
//    {
//        tVal( Ik, 0 ) = Ik;
//    }
//
//    //VecView ( tSolution, PETSC_VIEWER_STDOUT_WORLD);
//
//    VecGetValues ( tSolution, LHSValues.length(), tVal.data() , LHSValues.data() );
//
//    //VecDestroy( &tSolution );
//
//}
