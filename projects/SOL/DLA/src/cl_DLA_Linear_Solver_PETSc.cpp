/*
 * cl_DLA_Linear_Solver_PETSc.cpp
 *
 *  Created on: Dez 12, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Solver_PETSc.hpp"
#include "cl_Matrix_Vector_Factory.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

using namespace moris;
using namespace dla;

Linear_Solver_PETSc::Linear_Solver_PETSc()
{

    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::Linear_Solver_PETSc( Linear_Problem * aLinearSystem )
{
    mLinearSystem = aLinearSystem;

    //FIXME add rest
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::~Linear_Solver_PETSc()
{
    //KSPDestroy(&mPetscKSPProblem);
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_linear_problem( Linear_Problem * aLinearSystem )
{
    mLinearSystem = aLinearSystem;
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_parameters()
{
    // Create parameter list and set default values fo solver parameters

    // Set KSP type
    mParameterList.insert( "KSPType", std::string( KSPGMRES ) );

    // Set default preconditioner
    mParameterList.insert( "PCType", std::string( PCILU ) );

    // Sets maximal iters for KSP
    mParameterList.insert( "KSPMaxits", 1000 );

    // Sets KSP gmres restart
    mParameterList.insert( "KSPMGMRESRestart", 500 );

    // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES
    mParameterList.insert( "KSPGMRESHapTol", 1e-10 );

    // Sets tolerance for KSP
    mParameterList.insert( "KSPTol", 1e-10 );

    // Sets the number of levels of fill to use for ILU
    mParameterList.insert( "ILUFill", 0 );

    // Sets drop tolerance for ilu
    mParameterList.insert( "ILUTol", 1e-6 );

    // Set multigrid levels
    mParameterList.insert( "MultigridLevels", 3 );
}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system( )
{
    return 0;
}

//----------------------------------------------------------------------------------------
moris::sint Linear_Solver_PETSc::solve_linear_system(        Linear_Problem * aLinearSystem,
                                                      const moris::sint       aIter )
{

    // Create KSP and PC
    KSPCreate( PETSC_COMM_WORLD, &mPetscKSPProblem );
    KSPGetPC( mPetscKSPProblem, &mpc );

    this->set_solver_internal_parameters( );

    if ( ! strcmp(mParameterList.get< std::string >( "PCType" ).c_str(), "mg") )
    {
        this->build_multigrid_preconditioner( aLinearSystem );
    }

    KSPSetOperators( mPetscKSPProblem, aLinearSystem->get_matrix()->get_petsc_matrix(), aLinearSystem->get_matrix()->get_petsc_matrix() );
    KSPGMRESSetOrthogonalization( mPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );
    KSPSetFromOptions( mPetscKSPProblem );

//    aLinearSystem->get_free_solver_LHS()->read_vector_from_HDF5( "Exact_Sol_petsc.h5" );
//    aLinearSystem->get_free_solver_LHS()->print();

    aLinearSystem->get_solver_RHS()->save_vector_to_HDF5( "Res_vec.h5" );
//    aLinearSystem->get_solver_RHS()->print();

    // Solve System
    KSPSolve( mPetscKSPProblem, aLinearSystem->get_solver_RHS()->get_petsc_vector(), aLinearSystem->get_free_solver_LHS()->get_petsc_vector() );

    // Output
    //KSPView( mPetscKSPProblem, PETSC_VIEWER_STDOUT_WORLD );
    moris::sint Iter;
    KSPGetIterationNumber(mPetscKSPProblem, &Iter );
    std::cout<<Iter<<" Iterations"<<std::endl;

    return 0;
}

//----------------------------------------------------------------------------------------

void Linear_Solver_PETSc::build_multigrid_preconditioner( Linear_Problem * aLinearSystem )
{
    // Build multigrid operators
    aLinearSystem->get_solver_input()->build_multigrid_operators();

    // get multigrid operators
    moris::Cell< Sparse_Matrix * > tProlongationList = aLinearSystem->get_solver_input()->get_multigrid_operator_pointer()->get_prolongation_list();

    PetscInt tLevels = mParameterList.get< moris::sint >( "MultigridLevels" );

    PCMGSetLevels( mpc, tLevels, NULL );
    PCMGSetType( mpc, PC_MG_MULTIPLICATIVE );
    PCMGSetGalerkin( mpc, PETSC_TRUE );

    moris::Cell< Mat > tTransposeOperators( tProlongationList.size() );
    for ( moris::uint Ik = 0; Ik < tProlongationList.size(); Ik++ )
    {
        MatTranspose( tProlongationList( Ik )->get_petsc_matrix(), MAT_INITIAL_MATRIX, &tTransposeOperators( Ik ) );
    }

    moris::sint tCounter = 0;
    for ( moris::sint Ik = tLevels-1; Ik > 0; Ik-- )
    {
         PCMGSetInterpolation( mpc, Ik, tTransposeOperators( tCounter++ ) );
    }

    for ( moris::uint Ik = 0; Ik < tTransposeOperators.size(); Ik++ )
    {
        MatDestroy( &tTransposeOperators( Ik ) );
    }
//------------------------------------------
     KSP tPetscKSPCoarseSolve;
     PCMGGetCoarseSolve(mpc, &tPetscKSPCoarseSolve);
     KSPSetType(tPetscKSPCoarseSolve,KSPRICHARDSON);
     PetscInt maxitss=1000;
     KSPSetTolerances( tPetscKSPCoarseSolve, 1.e-10, PETSC_DEFAULT, PETSC_DEFAULT, maxitss );
     PC cpc;
     KSPGetPC(tPetscKSPCoarseSolve,&cpc);
     PCSetType(cpc,PCLU);
//----------------------------------------------------------------

     for (PetscInt k=1;k<tLevels;k++)
     {
         KSP dkspDown;
         PCMGGetSmootherDown(mpc,k,&dkspDown);
         PC dpcDown;
         KSPGetPC(dkspDown,&dpcDown);
         KSPSetType(dkspDown,KSPGMRES);                                                       // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
         moris::sint restart = 2;
         KSPGMRESSetRestart(dkspDown,restart);
         KSPSetTolerances(dkspDown,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,restart);         // NOTE maxitr=restart;
         PCSetType(dpcDown,PCJACOBI);                                                          // PCJACOBI, PCSOR for KSPCHEBYSHEV very good... Use KSPRICHARDSON for weighted Jacobi
     }

     for (PetscInt k=1;k<tLevels;k++)
     {
         KSP dkspUp;

         PCMGGetSmootherUp(mpc,k,&dkspUp);
         PC dpcUp;
         KSPGetPC(dkspUp,&dpcUp);
         KSPSetType(dkspUp,KSPGMRES);
//            KSPSetType(dkspUp,KSPGMRES);                                                                 // KSPCG, KSPGMRES, KSPCHEBYSHEV (VERY GOOD FOR SPD)
         moris::sint restart = 4;
         KSPGMRESSetRestart(dkspUp,restart);
         KSPSetTolerances(dkspUp,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,restart);                     // NOTE maxitr=restart;
         PCSetType(dpcUp,PCJACOBI);
     }
}

//----------------------------------------------------------------------------------------

void Linear_Solver_PETSc::set_solver_internal_parameters( )
{
        // Set KSP type
        KSPSetType( mPetscKSPProblem, mParameterList.get< std::string >( "KSPType" ).c_str() );
        KSPSetInitialGuessNonzero( mPetscKSPProblem, PETSC_TRUE );

        // Set maxits and tolerance for ksp
        KSPSetTolerances( mPetscKSPProblem, mParameterList.get< moris::real >( "KSPTol" ), PETSC_DEFAULT, PETSC_DEFAULT, mParameterList.get< moris::sint >( "KSPMaxits" ) );

        // Set Gmres restart
        KSPGMRESSetRestart( mPetscKSPProblem, mParameterList.get< moris::sint >( "KSPMGMRESRestart" ) );

        // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES.
        KSPGMRESSetHapTol( mPetscKSPProblem, mParameterList.get< moris::real >( "KSPGMRESHapTol" ) );

        // Set PC type
        PCSetType( mpc, mParameterList.get< std::string >( "PCType" ).c_str() );

        // Set levels of fill for ILU
        PCFactorSetLevels( mpc, mParameterList.get< moris::sint >( "ILUFill" ) );

        // Set drop tolerance for Ilu
        PCFactorSetDropTolerance( mpc, mParameterList.get< moris::real >( "ILUTol" ), PETSC_DEFAULT, PETSC_DEFAULT );

        PCSORSetOmega( mpc, 1 );

        PCSORSetIterations( mpc, 1 , 1 );
}


