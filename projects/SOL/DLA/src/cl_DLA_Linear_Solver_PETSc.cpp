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
#include "fn_PRM_SOL_Parameters.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

#include "cl_Tracer.hpp"

#ifdef MORIS_HAVE_SLEPC
#include "slepceps.h"
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

Linear_Solver_PETSc::Linear_Solver_PETSc()
{
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------

Linear_Solver_PETSc::Linear_Solver_PETSc( const moris::ParameterList aParameterlist )
        : Linear_Solver_Algorithm_Petsc( aParameterlist )
{
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::Linear_Solver_PETSc( Linear_Problem *aLinearSystem )
{
    mLinearSystem = aLinearSystem;

    // FIXME add rest
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------
Linear_Solver_PETSc::~Linear_Solver_PETSc()
{
    // KSPDestroy(&mPetscKSPProblem);
    //    PCDestroy(&mpc);
}

//----------------------------------------------------------------------------------------
void Linear_Solver_PETSc::set_solver_parameters()
{
    // Create parameter list and set default values for solver parameters

    // Set KSP type
    mParameterList.insert( "KSPType", std::string( "gmres" ) );

    // Set default preconditioner
    mParameterList.insert( "PCType", std::string( "ilu" ) );

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

    // Set multigrid levels
    mParameterList.insert( "ouput_eigenspectrum", (uint)0 );
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

    // Solve System
    KSPSolve(
            mPetscKSPProblem,
            static_cast< Vector_PETSc * >( aLinearSystem->get_solver_RHS() )->get_petsc_vector(),
            static_cast< Vector_PETSc * >( aLinearSystem->get_free_solver_LHS() )->get_petsc_vector() );

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
            NULL,
            0 );
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
    EPSSetOperators( eps, tBA, NULL );
    EPSSetProblemType( eps, EPS_HEP );
    EPSSetDimensions( eps, (PetscInt)tNumEigenValues, PETSC_DEFAULT, PETSC_DEFAULT );
    EPSSolve( eps );
    EPSGetTolerances( eps, &tTolerance, &tMaxIter );
#endif
}
