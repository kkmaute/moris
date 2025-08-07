/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_PETSc.cpp
 *
 */

#include "cl_DLA_Linear_Solver_PETSc.hpp"
#include "cl_DLA_Preconditioner_PETSc.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

#include "cl_Tracer.hpp"

#ifdef MORIS_HAVE_SLEPC
#include <slepceps.h>
#endif

using namespace moris;
using namespace dla;

//----------------------------------------------------------------------------------------

PetscErrorCode
fn_KSPMonitorResidual( KSP ksp, PetscInt n, PetscReal rnorm, void *dummy )
{
    MORIS_LOG_INFO( "KSP Iteration %d: Residual norm = %e", n, rnorm );
    return 0;
}

//----------------------------------------------------------------------------------------

Linear_Solver_PETSc::Linear_Solver_PETSc( const moris::Parameter_List& aParameterlist )
        : Linear_Solver_Algorithm_Petsc( aParameterlist )
{
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::Linear_Solver_PETSc( Linear_Problem *aLinearSystem )
        : Linear_Solver_Algorithm_Petsc( prm::create_linear_algorithm_parameter_list_petsc() )
{
    mLinearSystem = aLinearSystem;
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::~Linear_Solver_PETSc()
{
    // KSPDestroy(&mPetscKSPProblem);
    //    PCDestroy(&mpc);
}

//----------------------------------------------------------------------------------------
moris::sint
Linear_Solver_PETSc::solve_linear_system()
{
    MORIS_ERROR( false, "Linear_Solver_PETSc::solve_linear_system - function not implemented." );

    return 0;
}

//----------------------------------------------------------------------------------------
moris::sint
Linear_Solver_PETSc::solve_linear_system(
        Linear_Problem   *aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "PETSc", "Solve" );

    // Create KSP
    KSPCreate( PETSC_COMM_WORLD, &mPetscKSPProblem );

    // Set matrices for linear system and for preconditioner
    KSPSetOperators( mPetscKSPProblem,
            aLinearSystem->get_matrix()->get_petsc_matrix(),
            aLinearSystem->get_matrix()->get_petsc_matrix() );

    // set solver interface (used by preconditioners)
    mSolverInterface = aLinearSystem->get_solver_input();

    // construct solver and preconditioner
    this->construct_solver_and_preconditioner( aLinearSystem );

    // for debugging: print matrix, rhs, and lhs
    // MatView( aLinearSystem->get_matrix()->get_petsc_matrix(), PETSC_VIEWER_STDOUT_WORLD );
    // VecView( static_cast< Vector_PETSc * >( aLinearSystem->get_solver_RHS() )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );
    // VecView( static_cast< Vector_PETSc * >( aLinearSystem->get_free_solver_LHS() )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );

    this->compute_eigenspectrum( aLinearSystem );

    Mat tRHSVecs = static_cast< MultiVector_PETSc* >( aLinearSystem->get_solver_RHS() )->get_petsc_vector();
    Mat tLHSVecs = static_cast< MultiVector_PETSc* >( aLinearSystem->get_free_solver_LHS() )->get_petsc_vector();

    for ( uint iNumRHS = 0; iNumRHS < mSolverInterface->get_num_rhs(); iNumRHS++ )
    {
        Vec tRHSVec, tLHSVec;
        MatDenseGetColumnVec( tRHSVecs, iNumRHS, &tRHSVec );
        MatDenseGetColumnVec( tLHSVecs, iNumRHS, &tLHSVec );

        VecAssemblyBegin( tRHSVec );
        VecAssemblyEnd( tRHSVec );
        VecAssemblyBegin( tLHSVec );
        VecAssemblyEnd( tLHSVec );

        KSPSolve( mPetscKSPProblem,tRHSVec,tLHSVec );
        MatDenseRestoreColumnVec( tRHSVecs, iNumRHS, &tRHSVec );
        MatDenseRestoreColumnVec( tLHSVecs, iNumRHS, &tLHSVec );
    }

    // for debugging: print lhs after solve
    // VecView( static_cast< Vector_PETSc * >( aLinearSystem->get_free_solver_LHS() )->get_petsc_vector(), PETSC_VIEWER_STDOUT_WORLD );

    mSolverInterface = nullptr;

    KSPDestroy( &mPetscKSPProblem );

    return 0;
}

//----------------------------------------------------------------------------------------

void Linear_Solver_PETSc::set_solver_analysis_options()
{
    KSPMonitorSet( mPetscKSPProblem,
            fn_KSPMonitorResidual,
            nullptr,
            nullptr );
}

//----------------------------------------------------------------------------------------

void Linear_Solver_PETSc::construct_solver_and_preconditioner( Linear_Problem *aLinearSystem )
{

    if ( !strcmp( mParameterList.get< std::string >( "KSPType" ).c_str(), "preonly" ) )
    {
        // write solver to log file
        MORIS_LOG_INFO( "KSP Solver: preonly" );

        // set solver
        KSPSetType( mPetscKSPProblem, KSPPREONLY );
        mPreconditioner->build_preconditioner( aLinearSystem, mPetscKSPProblem );    //   KSPGetPC( mPetscKSPProblem, mPreconditioner->get_pc() );
    }
    // set iterative solver: kspgmres
    else if (                                                                                //
            !strcmp( mParameterList.get< std::string >( "KSPType" ).c_str(), "gmres" ) ||    //
            !strcmp( mParameterList.get< std::string >( "KSPType" ).c_str(), "fgmres" ) )
    {

        // write solver to log file
        MORIS_LOG_INFO( "KSP Solver: %s", mParameterList.get< std::string >( "KSPType" ).c_str() );

        // set solver type
        KSPSetType( mPetscKSPProblem, mParameterList.get< std::string >( "KSPType" ).c_str() );

        // use initial guess
        // KSPSetInitialGuessNonzero( mPetscKSPProblem, PETSC_TRUE );

        // set orthogonalization method for gmres
        KSPGMRESSetOrthogonalization( mPetscKSPProblem, KSPGMRESModifiedGramSchmidtOrthogonalization );

        // Set maxits and tolerance for ksp
        KSPSetTolerances(
                mPetscKSPProblem,
                mParameterList.get< moris::real >( "KSPTol" ),
                PETSC_DEFAULT,
                PETSC_DEFAULT,
                mParameterList.get< moris::sint >( "KSPMaxits" ) );

        // Set Gmres restart
        KSPGMRESSetRestart( mPetscKSPProblem, mParameterList.get< moris::sint >( "KSPMGMRESRestart" ) );

        // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES.
        KSPGMRESSetHapTol( mPetscKSPProblem, mParameterList.get< moris::real >( "KSPGMRESHapTol" ) );

        // initialize preconditioner
        if ( mPreconditioner != nullptr )
        {
            mPreconditioner->build_preconditioner( aLinearSystem, mPetscKSPProblem );
        }
    }
    else if ( !strcmp( mParameterList.get< std::string >( "KSPType" ).c_str(), "stcg" ) )
    {
        // write solver to log file
        MORIS_LOG_INFO( "KSP Solver: stcg - for use in trust region solver only" );

        // set solver type
        KSPSetType( mPetscKSPProblem, mParameterList.get< std::string >( "KSPType" ).c_str() );

        // use initial guess
        // KSPSetInitialGuessNonzero( mPetscKSPProblem, PETSC_TRUE );

        // Set maxits and tolerance for ksp
        KSPSetTolerances(
                mPetscKSPProblem,
                mParameterList.get< moris::real >( "KSPTol" ),
                PETSC_DEFAULT,
                PETSC_DEFAULT,
                mParameterList.get< moris::sint >( "KSPMaxits" ) );

        // Set the trust region size
        KSPCGSetRadius( mPetscKSPProblem, mTrSize );

        // initialize preconditioner
        if ( mPreconditioner != nullptr )
        {
           mPreconditioner->build_preconditioner( aLinearSystem, mPetscKSPProblem );
        }
    }
    else
    {
        MORIS_ERROR( false, "Linear_Solver_PETSc::construct_solver_and_preconditioner - unsupported KSP type." );
    }

    // set convergence options
    this->set_solver_analysis_options();

    // finalize solver setup
    KSPSetFromOptions( mPetscKSPProblem );

    // finalize solver setup
    KSPSetUp( mPetscKSPProblem );

    // KSPSetComputeEigenvalues(mPetscKSPProblem,PETSC_TRUE);

    // for debugging: print solver setup
    KSPView( mPetscKSPProblem, PETSC_VIEWER_STDOUT_WORLD );
}

void Linear_Solver_PETSc::compute_eigenspectrum( Linear_Problem *aLinearSystem )
{
    uint tNumEigenValues = mParameterList.get< uint >( "ouput_eigenspectrum" );
    if ( tNumEigenValues == 0 )
    {
        return;
    }
#ifdef MORIS_HAVE_SLEPC

    // declare the explict precondioied matrix
    ::Mat tBA;

    // get the matrix type fo A
    Mat tMat = aLinearSystem->get_matrix()->get_petsc_matrix();

    // get the matrix type
    MatType tMatType;
    MatGetType( tMat, &tMatType );

    // compute the explicit operator
    KSPComputeOperator( mPetscKSPProblem, tMatType, &tBA );

    // set up the eigen problem
    EPS         eps;
    moris::real tTolerance;
    moris::sint tMaxIter;
    SlepcInitializeNoArguments();
    EPSCreate( PETSC_COMM_WORLD, &eps );
    EPSSetOperators( eps, tBA, nullptr );
    EPSSetProblemType( eps, EPS_HEP );
    EPSSetDimensions( eps, (PetscInt)tNumEigenValues, PETSC_DEFAULT, PETSC_DEFAULT );
    EPSSolve( eps );
    EPSGetTolerances( eps, &tTolerance, &tMaxIter );
#endif
}
