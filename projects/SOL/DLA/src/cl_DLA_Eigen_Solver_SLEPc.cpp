/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Eigen_Solver_SLEPc.cpp
 *
 */

#include "cl_DLA_Eigen_Solver_SLEPc.hpp"
#include "cl_DLA_Preconditioner_PETSc.hpp"
#include "cl_SOL_Matrix_Vector_Factory.hpp"
#include "fn_PRM_SOL_Parameters.hpp"

#include <petscksp.h>
#include <petscdm.h>
#include <petscdmda.h>
#include "petscmat.h"

#include <string>

#include "cl_Tracer.hpp"
#include "slepceps.h"

#include "cl_Vector_PETSc_Multi.hpp"

using namespace moris;
using namespace dla;

//----------------------------------------------------------------------------------------

PetscErrorCode
fn_EPSMonitorResidual( EPS eps, PetscInt its, PetscInt nconv, PetscScalar *eigr, PetscScalar *eigi, PetscReal *errest, PetscInt nest, void *ctx )
{
    MORIS_LOG_INFO( "EPS Iteration %d: number of converged eigenpairs = %d", its, nconv );

    // collect all values in a string
    std::ostringstream tOss;
    for ( PetscInt i = 0; i < nconv; ++i )
    {
        tOss << errest[ i ] << '\t';
    }

    // output the string
    MORIS_LOG_INFO( "%s", tOss.str().c_str() );
    return 0;
}

//----------------------------------------------------------------------------------------

// Eigen_Solver_SLEPc::Eigen_Solver_SLEPc()
// {
//      SlepcInitializeNoArguments();
//     mParameterList = prm::create_slepc_algorithm_parameter_list();
// }

//----------------------------------------------------------------------------------------

Eigen_Solver_SLEPc::Eigen_Solver_SLEPc( const moris::Parameter_List aParameterlist )
        : Linear_Solver_Algorithm_Petsc( aParameterlist )
{
    SlepcInitializeNoArguments();

    mStringToEPSWhich = {
        { "LM", EPS_LARGEST_MAGNITUDE },
        { "SM", EPS_SMALLEST_MAGNITUDE },
        { "LR", EPS_LARGEST_REAL },
        { "SR", EPS_SMALLEST_REAL },
        { "LI", EPS_LARGEST_IMAGINARY },
        { "SI", EPS_SMALLEST_IMAGINARY },
        { "TM", EPS_TARGET_MAGNITUDE },
        { "TR", EPS_TARGET_REAL },
        { "TI", EPS_TARGET_IMAGINARY },
        { "A", EPS_ALL },
        { "U", EPS_WHICH_USER }
    };

    mMorisSTToSTType = {
        { "shell", STSHELL },
        { "precond", STPRECOND },
        { "filter", STFILTER },
        { "shift", STSHIFT },
        { "shift_invert", STSINVERT },
        { "cayley", STCAYLEY }
    };

    // map for the power method shift type
    mMorisToPowerType = {
        { "constant", EPS_POWER_SHIFT_CONSTANT },
        { "rayleigh", EPS_POWER_SHIFT_RAYLEIGH },
        { "spectral", EPS_POWER_SHIFT_WILKINSON }
    };

    mConvergenceReasonToString = {
        { EPS_CONVERGED_TOL, "EPS_CONVERGED_TOL" },                      // converged up to tolerance
        { EPS_CONVERGED_ITERATING, "EPS_CONVERGED_ITERATING" },          // 0 value
        { EPS_CONVERGED_USER, "EPS_CONVERGED_USER" },                    // converged due to a user-defined condition
        { EPS_DIVERGED_BREAKDOWN, "EPS_DIVERGED_BREAKDOWN" },            // generic breakdown in method
        { EPS_DIVERGED_SYMMETRY_LOST, "EPS_DIVERGED_SYMMETRY_LOST" },    // pseudo-Lanczos was not able to keep symmetry
        { EPS_DIVERGED_ITS, "EPS_DIVERGED_ITS" }                         // required more than max_it iterations to reach convergence
    };
}

//----------------------------------------------------------------------------------------
Eigen_Solver_SLEPc::~Eigen_Solver_SLEPc()
{
    // KSPDestroy(&mPetscKSPProblem);
    //    PCDestroy(&mpc);
}


//----------------------------------------------------------------------------------------
moris::sint
Eigen_Solver_SLEPc::solve_linear_system()
{
    MORIS_ERROR( false, "Eigen_Solver_SLEPc::solve_linear_system - function not implemented." );
    return 0;
}

//----------------------------------------------------------------------------------------
moris::sint
Eigen_Solver_SLEPc::solve_linear_system(
        Linear_Problem   *aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "SLEPc", "Solve" );


    // cteate the eigen problem solve object
    EPSCreate( PETSC_COMM_WORLD, &mEps );
    EPSProblemType tProblemType = this->determine_problem_type( aLinearSystem );
    EPSSetProblemType( mEps, tProblemType );

    // detrmine what solver to use , the defualt is power
    this->set_eps_type_and_params();

    // set whcih eigenvalues to solve for
    EPSSetWhichEigenpairs( mEps, mStringToEPSWhich.find( mParameterList.get< std::string >( "Which" ) )->second );

    // request number of eigenvalues from parameterlist
    moris::sint tNumEigVals = mParameterList.get< moris::sint >( "Num_Eig_Vals" );
    EPSSetDimensions( mEps, (PetscInt)tNumEigVals, PETSC_DEFAULT, PETSC_DEFAULT );

    if ( mParameterList.get< bool >( "Verbosity" ) ) EPSMonitorSet( mEps, fn_EPSMonitorResidual, NULL, 0 );

    EPSGetST( mEps, &mSt );
    if ( mParameterList.get< std::string >( "STType" ) == "shift_invert" )
    {
        EPSSetTarget( mEps, mParameterList.get< real >( "ShiftValue" ) );
        EPSSetWhichEigenpairs( mEps, EPS_TARGET_MAGNITUDE );
    }
    else
    {
        STSetShift( mSt, mParameterList.get< real >( "ShiftValue" ) );
    }
    STSetType( mSt, mMorisSTToSTType.find( mParameterList.get< std::string >( "STType" ) )->second );


    // set the sub linear solver and preconditioner options
    this->set_sublinear_solver_and_preconditioner( aLinearSystem );

    real tTolerance = mParameterList.get< real >( "Convergence_Tolerance" );
    uint tMaxIter   = mParameterList.get< uint >( "max_iter" );
    EPSSetTolerances( mEps, tTolerance, tMaxIter );


    EPSSolve( mEps );

    // Get the convergence reason
    EPSGetConvergedReason( mEps, &reason );
    MORIS_LOG_INFO( " Finished - converged reason = %s \n", mConvergenceReasonToString.find( reason )->second.c_str() );

    EPSView( mEps, PETSC_VIEWER_STDOUT_WORLD );
    if ( mParameterList.get< bool >( "Verbosity" ) ) print_slepc_determined_solver_paramaters();

    moris::moris_id tNumConvergedEigVals;
    EPSGetConverged( mEps, &tNumConvergedEigVals );
    MORIS_LOG_INFO( "Number of converged eigenvalues: %d", tNumConvergedEigVals );

    // iterate over the number of converged eigenvalues and get the corresponding eigenvalues
    for ( moris_id iEigenIndex = 0; iEigenIndex < tNumConvergedEigVals; iEigenIndex++ )
    {
        real tEigenValueReal, tEigenValueImag, tError;
        EPSGetEigenvalue( mEps, iEigenIndex, &tEigenValueReal, &tEigenValueImag );
        EPSComputeError( mEps, iEigenIndex, EPS_ERROR_RELATIVE, &tError );

        MORIS_LOG_INFO( "Eigenvalue %d : %f + %fi , Error : %f", iEigenIndex, tEigenValueReal, tEigenValueImag, tError );

        mEigenValues.push_back( tEigenValueReal );

        if ( aLinearSystem->get_solver_input() not_eq nullptr )
        {
            std::shared_ptr< Vector< real > > &tEigenValues = aLinearSystem->get_solver_input()->get_eigen_values();
            tEigenValues->push_back( tEigenValueReal );
        }
    }

    if ( !mParameterList.get< bool >( "Update_Flag" ) ) return 0;

    Vec tSourceVec;    // petsc vector
    for ( moris_id iEigenIndex = 0; iEigenIndex < tNumEigVals; iEigenIndex++ )
    {
        sol::Dist_Vector  *tDistVec           = aLinearSystem->get_solver_input()->get_eigen_solution_vector();
        MultiVector_PETSc *tDestinationVector = static_cast< MultiVector_PETSc * >( tDistVec );

        // get the eigen vector
        MatCreateVecs( aLinearSystem->get_matrix()->get_petsc_matrix(), NULL, &tSourceVec );
        EPSGetEigenvector( mEps, iEigenIndex, tSourceVec, NULL );

        tDestinationVector->import_local_to_global( tSourceVec, iEigenIndex );
    }


    return 0;
}

//----------------------------------------------------------------------------------------

EPSProblemType
Eigen_Solver_SLEPc::determine_problem_type( Linear_Problem *aLinearSystem )
{
    // get the matrix type
    // determine if the right hand side is the mass matrix or identity matrix
    std::string tRHSType = aLinearSystem->get_rhs_matrix_type();

    // determine if the problem is symmetric
    bool tAssumeSymmetric = mParameterList.get< bool >( "is_symmetric" );

    // determine the problem type based on symmetry and right hand side type
    if ( tRHSType == "IdentityMat" )
    {
        EPSSetOperators( mEps, aLinearSystem->get_matrix()->get_petsc_matrix(), NULL );
        return tAssumeSymmetric ? EPS_HEP : EPS_NHEP;
    }
    else if ( tRHSType == "MassMat" )
    {
        aLinearSystem->assemble_rhs_matrix();
        EPSSetOperators( mEps, aLinearSystem->get_matrix()->get_petsc_matrix(), aLinearSystem->get_mass_matrix()->get_petsc_matrix() );
        return tAssumeSymmetric ? EPS_GHEP : EPS_GNHEP;
    }
    else if ( tRHSType == "GeomStiffMat" )
    {
        aLinearSystem->assemble_rhs_matrix();
        EPSSetOperators( mEps, aLinearSystem->get_matrix()->get_petsc_matrix(), aLinearSystem->get_mass_matrix()->get_petsc_matrix() );
        return tAssumeSymmetric ? EPS_GHEP : EPS_GNHEP;
    }
    else
    {
        // Handle other cases here
        MORIS_ERROR( false, "Eigen_Solver_SLEPc::determine_problem_type - problem type not recognized." );
        return EPS_HEP;
    }
}

//----------------------------------------------------------------------------------------

extern PetscErrorCode
fn_KSPMonitorResidual( KSP ksp, PetscInt n, PetscReal rnorm, void *dummy );

void Eigen_Solver_SLEPc::set_sublinear_solver_options( const Parameter_List *aParameterlistsubSolver, const Parameter_List *aParameterlistPreconditioner )
{
    // petsc context thus cast into petsc linear solver algorthim to extract the petsc object
    mSubSolverParameterlist               = aParameterlistsubSolver;
    mSubSolverPreconditionerParameterlist = aParameterlistPreconditioner;
}

void Eigen_Solver_SLEPc::set_sublinear_solver_and_preconditioner( Linear_Problem *aLinearSystem )
{
    if ( mSubSolverParameterlist->size() == 0 )
    {
        return;
    }
    STGetKSP( mSt, &mKsp );
    MORIS_LOG_INFO( "KSP Solver: %s", mSubSolverParameterlist->get< std::string >( "KSPType" ).c_str() );

    // set the KSP options
    // set the options as necessary , need to check for the precondioner
    if ( !strcmp( mSubSolverParameterlist->get< std::string >( "KSPType" ).c_str(), "preonly" ) )
    {

        // set solver
        KSPSetType( mKsp, KSPPREONLY );
    }
    // set iterative solver: kspgmres
    else if (                                                                                          //
            !strcmp( mSubSolverParameterlist->get< std::string >( "KSPType" ).c_str(), "gmres" ) ||    //
            !strcmp( mSubSolverParameterlist->get< std::string >( "KSPType" ).c_str(), "fgmres" ) )
    {

        // set solver type
        KSPSetType( mKsp, mSubSolverParameterlist->get< std::string >( "KSPType" ).c_str() );

        // use initial guess
        // KSPSetInitialGuessNonzero( mKsp, PETSC_TRUE );

        // set orthogonalization method for gmres
        KSPGMRESSetOrthogonalization( mKsp, KSPGMRESModifiedGramSchmidtOrthogonalization );

        // Set maxits and tolerance for ksp
        KSPSetTolerances(
                mKsp,
                mSubSolverParameterlist->get< moris::real >( "KSPTol" ),
                PETSC_DEFAULT,
                PETSC_DEFAULT,
                mSubSolverParameterlist->get< moris::sint >( "KSPMaxits" ) );

        // Set Gmres restart
        KSPGMRESSetRestart( mKsp, mSubSolverParameterlist->get< moris::sint >( "KSPMGMRESRestart" ) );

        // Sets tolerance for determining happy breakdown in GMRES, FGMRES and LGMRES.
        KSPGMRESSetHapTol( mKsp, mSubSolverParameterlist->get< moris::real >( "KSPGMRESHapTol" ) );
    }

    if ( mParameterList.get< bool >( "Verbosity" ) ) KSPMonitorSet( mKsp, fn_KSPMonitorResidual, NULL, 0 );

    // set the preconditioner
    this->build_preconditioner( aLinearSystem );

    // for debugging: print solver setup
    if ( mParameterList.get< bool >( "Verbosity" ) ) KSPView( mKsp, PETSC_VIEWER_STDOUT_WORLD );
}

//----------------------------------------------------------------------------------------
void Eigen_Solver_SLEPc::build_preconditioner( Linear_Problem *aLinearProblem )
{
    // get preconditioner
    KSPGetPC( mKsp, &mPc );

    // set direct solver: superlu-dist
    if ( !strcmp( mSubSolverPreconditionerParameterlist->get< std::string >( "PCType" ).c_str(), "superlu-dist" ) )
    {
        // write solver to log file
        MORIS_LOG_INFO( "PC Solver: superlu-dist" );

        // set LU preconditioner
        PCSetType( mPc, PCLU );

        // set factorization method in preconditioner
        PCFactorSetMatSolverType( mPc, MATSOLVERSUPERLU_DIST );

        STGetOperator( mSt, NULL );

        // set up the package to call for the factorization
        PCFactorSetUpMatSolverType( mPc );
    }

    // set direct solver: mumps
    else if ( !strcmp( mSubSolverPreconditionerParameterlist->get< std::string >( "PCType" ).c_str(), "mumps" ) )
    {
#ifdef MORIS_USE_MUMPS
        // write solver to log file
        MORIS_LOG_INFO( "PC Solver: mumps" );

        // set LU preconditioner assuming that system is non-symmetric
        PCSetType( mPc, PCLU );

        // set factorization method in preconditioner
        PCFactorSetMatSolverType( mPc, MATSOLVERMUMPS );

        STGetOperator( mSt, NULL );

        // set up the package to call for the factorization
        PCFactorSetUpMatSolverType( mPc );

        // get the factored matrix F from the preconditioner context
        Mat F;
        PCFactorGetMatrix( mPc, &F );

        // set MUMPS integer control parameters ICNTL to be passed to
        // MUMPS.  Setting entry 7 of MUMPS ICNTL array (of size 40) to a value
        // of 2. This sets use of Approximate Minimum Fill (AMF)
        PetscInt ival = 2, icntl = 7;

        // pass control parameters to MUMPS

        MatMumpsSetIcntl( F, icntl, ival );
#else
        MORIS_ERROR( false,
                "Linear_Solver_PETSc::construct_solver_and_preconditioner - MORIS installed without support for MUMPS" );
#endif
    }

    else if ( !strcmp( mSubSolverPreconditionerParameterlist->get< std::string >( "PCType" ).c_str(), "none" ) )
    {
        // Set PC type to none
        PCSetType( mPc, "none" );
    }
    else
    {
        MORIS_ERROR( false,
                "Linear_Solver_PETSc::construct_solver_and_preconditioner - no valid preconditioner was found." );
    }
}

//----------------------------------------------------------------------------------------
void Eigen_Solver_SLEPc::set_eps_type_and_params()
{
    std::string tSolverType = mParameterList.get< std::string >( "Eigen_Algorithm" );
    bool        tUseDefualt = mParameterList.get< bool >( "use_slepc_defualt_params" );

    if ( tSolverType == "power" )
    {
        EPSSetType( mEps, tSolverType.c_str() );
        if ( tUseDefualt ) return;

        EPSPowerSetShiftType( mEps, mMorisToPowerType.find( mParameterList.get< std::string >( "shift_type" ) )->second );
    }
    else if ( tSolverType == "gd" )
    {
        // refer to fn_PRM_SOL_Parameters.hpp for the default values
        EPSSetType( mEps, tSolverType.c_str() );
        if ( tUseDefualt ) return;

        EPSGDSetKrylovStart( mEps, PetscBool( mParameterList.get< bool >( "krylov_start" ) ) );
        EPSGDSetBlockSize( mEps, mParameterList.get< uint >( "block_size" ) );
        EPSGDSetRestart( mEps,    //
                mParameterList.get< uint >( "number_of_vectors_initial_search_subspace" ),
                mParameterList.get< uint >( "number_of_vectors_after_restart" ) );
        EPSGDSetInitialSize( mEps,    //
                mParameterList.get< sint >( "number_of_vectors_saved_from_previous_restart" ) );
        EPSGDSetBOrth( mEps, PetscBool( mParameterList.get< bool >( "use_B_ortho" ) ) );
        EPSGDSetDoubleExpansion( mEps, PetscBool( mParameterList.get< bool >( "dynamic" ) ) );
    }
    else if ( tSolverType == "krylovschur" )
    {
        EPSSetType( mEps, tSolverType.c_str() );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "arnoldi" )
    {
        EPSSetType( mEps, EPSARNOLDI );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "lapack" )
    {
        EPSSetType( mEps, EPSLAPACK );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "jd" )
    {
        EPSSetType( mEps, EPSJD );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "lapack" )
    {
        EPSSetType( mEps, EPSLAPACK );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "arpack" )
    {
        EPSSetType( mEps, EPSARPACK );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "lapack" )
    {
        EPSSetType( mEps, EPSLAPACK );
        if ( tUseDefualt ) return;
    }
    else if ( tSolverType == "lyapii" )
    {
        EPSSetType( mEps, EPSLYAPII );
        if ( tUseDefualt ) return;
    }
    else
    {
        MORIS_ERROR( false, "Eigen_Solver_SLEPc::set_eps_type_and_params - solver type not recognized." );
    }

    MORIS_LOG_WARNING( "Eigen_Solver_SLEPc::set_eps_type_and_params - slepc defualt parametrs are being chosen" );
}

void Eigen_Solver_SLEPc::print_slepc_determined_solver_paramaters()
{
    std::string tSolverType = mParameterList.get< std::string >( "Eigen_Algorithm" );

    if ( tSolverType == "power" )
    {
        EPSPowerShiftType tShiftType;
        EPSPowerGetShiftType( mEps, &tShiftType );    // mMorisToPowerType.find( mParameterList.get<std::string>("shift_type") )->second );


        std::string shiftTypeString;
        for ( const auto &pair : mMorisToPowerType )
        {
            if ( pair.second == tShiftType )
            {
                shiftTypeString = pair.first;
                break;
            }
        }
        MORIS_LOG_INFO( "Shift type: %s", shiftTypeString.c_str() );
    }
    else if ( tSolverType == "gd" )
    {
        PetscBool op, orth;
        PetscInt  opi, opi0;
        EPSGDGetKrylovStart( mEps, &op );
        MORIS_LOG_INFO( "Krylov start: %d\n", op );
        EPSGDGetBOrth( mEps, &orth );
        MORIS_LOG_INFO( "BOrth: %d\n", orth );
        EPSGDGetBlockSize( mEps, &opi );
        MORIS_LOG_INFO( "Block size: %d\n", opi );
        EPSGDGetRestart( mEps, &opi, &opi0 );
        MORIS_LOG_INFO( "Restart: %d %d\n", opi, opi0 );
        EPSGDGetInitialSize( mEps, &opi );
        MORIS_LOG_INFO( "Initial size: %d\n", opi );
        EPSGDGetDoubleExpansion( mEps, &op );
        MORIS_LOG_INFO( "Double expansion(dynamic): %d\n", op );
    }
    else if ( tSolverType == "krylovschur" )
    {
        MORIS_LOG_INFO( "Eiegn solver defualt paramater is not implemenetd use slepc documenation!" );
        EPSView( mEps, PETSC_VIEWER_STDOUT_WORLD );
    }
    else
    {
        MORIS_LOG_INFO( "Eiegn solver defualt paramater is not implemenetd use slepc documenation!" );
        EPSView( mEps, PETSC_VIEWER_STDOUT_WORLD );
    }
}

//----------------------------------------------------------------------------------------

Vector< real > const &
Eigen_Solver_SLEPc::get_eigenvalues() const
{
    return mEigenValues;
}
