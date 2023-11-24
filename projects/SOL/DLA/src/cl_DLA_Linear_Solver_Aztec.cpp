/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Aztec.cpp
 *
 */
#include "cl_DLA_Linear_Solver_Aztec.hpp"

// TPL header files
#include "Epetra_ConfigDefs.h"
#include "AztecOO_ConfigDefs.h"
#include "AztecOO_ConditionNumber.h"

#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#include "cl_DLA_Linear_Problem.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

// Teuchos
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

// detailed logging
#include "cl_Tracer.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_Aztec::Linear_Solver_Aztec()
{
    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------

Linear_Solver_Aztec::Linear_Solver_Aztec( const moris::ParameterList aParameterlist )
        : Linear_Solver_Algorithm_Trilinos( aParameterlist )
{
}

//----------------------------------------------------------------------------------------

Linear_Solver_Aztec::Linear_Solver_Aztec( Linear_Problem* aLinearSystem )
{
    // store linear problem for building external preconditioner
    mLinearSystem = aLinearSystem;

    // Set matrix. solution vector and RHS
    mEpetraProblem.SetOperator( mLinearSystem->get_matrix()->get_matrix() );

    mEpetraProblem.SetRHS( static_cast< Vector_Epetra* >(
            mLinearSystem->get_solver_RHS() )
                                   ->get_epetra_vector() );

    mEpetraProblem.SetLHS( static_cast< Vector_Epetra* >(
            mLinearSystem->get_free_solver_LHS() )
                                   ->get_epetra_vector() );

    this->set_solver_parameters();
}

//----------------------------------------------------------------------------------------

Linear_Solver_Aztec::~Linear_Solver_Aztec()
{
    delete mAztecSolver;
    mAztecSolver = nullptr;
}

//----------------------------------------------------------------------------------------

void
Linear_Solver_Aztec::set_solver_parameters()
{
    mParameterList = prm::create_linear_algorithm_parameter_list_aztec();
}

//----------------------------------------------------------------------------------------

bool
Linear_Solver_Aztec::build_external_preconditioner( const sint& aIter )
{
    // if it is not null
    if ( mPreconditioner )
    {
        // build preconditioner
        mPreconditioner->build( mLinearSystem, aIter );

        // check if preconditioner exists
        return mPreconditioner->exists();
    }

    return false;
}

// -----------------------------------------------------------------------------------

moris::sint
Linear_Solver_Aztec::solve_linear_system()
{
    Tracer tTracer( "LinearAlgorithm", "Aztec", "Solve" );

    // Construct solver with linear system
    mAztecSolver = new AztecOO( mEpetraProblem );

    // Check that problem has only one RHS
    Epetra_MultiVector* tLHS = mAztecSolver->GetLHS();

    MORIS_ERROR( tLHS->NumVectors() == 1,
            "AZTEC interface cannot be used for multiple RHS.\n" );

    moris::sint error = 0;

    // Set all Aztec options
    this->set_solver_internal_parameters();

    moris::sint tMaxIt = mParameterList.get< moris::sint >( "AZ_max_iter" );

    moris::real tRelRes = mParameterList.get< moris::real >( "rel_residual" );

    // Build preconditioner based on input parameters
    bool tPreconditionerIsBuilt = this->build_external_preconditioner();

    // Set external preconditioner if exists
    if ( tPreconditionerIsBuilt )
    {
        mAztecSolver->SetPrecOperator( mPreconditioner->get_operator().get() );
    }

    // Solve the linear system
    error = mAztecSolver->Iterate( tMaxIt, tRelRes );

    // Get linear solution info
    mSolNumIters       = mAztecSolver->NumIters();
    mSolTrueResidual   = mAztecSolver->TrueResidual();
    mSolScaledResidual = mAztecSolver->ScaledResidual();
    mSolTime           = mAztecSolver->SolveTime();

    // log linear solver iterations
    MORIS_LOG_SPEC( "LinearSolverIterations", mSolNumIters );

    // compute exact residuals
    MORIS_LOG_SPEC( "LinearResidualNorm", mSolTrueResidual );

    delete mAztecSolver;
    mAztecSolver = nullptr;

    return error;
}

//----------------------------------------------------------------------------------------

moris::sint
Linear_Solver_Aztec::solve_linear_system(
        Linear_Problem*   aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "Aztec", "Solve" );

    // set linear system
    mLinearSystem = aLinearSystem;

    // Set matrix in linear system
    mEpetraProblem.SetOperator( mLinearSystem->get_matrix()->get_matrix() );

    // Get LHS and RHS vectors
    sol::Dist_Vector* tRHS = mLinearSystem->get_solver_RHS();
    sol::Dist_Vector* tLHS = mLinearSystem->get_free_solver_LHS();

    // Determine the number of RHS and LHS
    uint tNumRHS = tRHS->get_num_vectors();
    uint tNumLHS = tLHS->get_num_vectors();

    MORIS_ERROR( tNumRHS == tNumLHS, "Number of LHS does not match number of RHS" );

    // Get underlying Eptra vectors of RHS and LHS
    Epetra_MultiVector* tRHSepetra = static_cast< Vector_Epetra* >( tRHS )->get_epetra_vector();
    Epetra_MultiVector* tLHSepetra = static_cast< Vector_Epetra* >( tLHS )->get_epetra_vector();

    // get basic solver parameters
    moris::sint tMaxIt  = mParameterList.get< moris::sint >( "AZ_max_iter" );
    moris::real tRelRes = mParameterList.get< moris::real >( "rel_residual" );

    // initialize and build external preconditioner based on input parameters
    bool tPreconditionerIsBuilt =  this->build_external_preconditioner( aIter );

    // initialize error flag
    moris::sint error = 0;

    // Loop over all RHS
    for ( uint ir = 0; ir < tNumRHS; ++ir )
    {
        // Get vectors of current RHSunderlying eptra vectors of RHS and LHS
        Epetra_MultiVector* tCurrentRHSepetra = ( *tRHSepetra )( ir );
        Epetra_MultiVector* tCurrentLHSepetra = ( *tLHSepetra )( ir );

        // Build linear problem and create solver for first RHS; otherwise just reset LHS and RHS
        if ( ir == 0 )
        {
            // Set matrix. solution vector and RHS
            mEpetraProblem.SetRHS( tCurrentRHSepetra );
            mEpetraProblem.SetLHS( tCurrentLHSepetra );

            mAztecSolver = new AztecOO( mEpetraProblem );

            // Set all Aztec options based on default and user input
            this->set_solver_internal_parameters();

            // Overwrite options for multiple RHS
            mAztecSolver->SetAztecOption( AZ_pre_calc, AZ_calc );
            if ( tNumRHS > 1 )
            {
                mAztecSolver->SetAztecOption( AZ_keep_info, 1 );
            }

            // Set external preconditioner if exists
            if ( tPreconditionerIsBuilt )
            {
                mAztecSolver->SetPrecOperator( mPreconditioner->get_operator().get() );
            }
        }
        else
        {
            mAztecSolver->SetRHS( tCurrentRHSepetra );
            mAztecSolver->SetLHS( tCurrentLHSepetra );

            // Overwrite options for multiple RHS
            mAztecSolver->SetAztecOption( AZ_pre_calc, AZ_reuse );
        }

        // Solve the linear system
        error += mAztecSolver->Iterate( tMaxIt, tRelRes );

        // Get linear solution info
        mSolNumIters       = mAztecSolver->NumIters();
        mSolTrueResidual   = mAztecSolver->TrueResidual();
        mSolScaledResidual = mAztecSolver->ScaledResidual();
        mSolTime           = mAztecSolver->SolveTime();

        const double* tStatus = mAztecSolver->GetAztecStatus();
        MORIS_LOG_SPEC( "Condition Number of Operator: ", tStatus[ AZ_condnum ] );
    }

    // log linear solver iterations
    MORIS_LOG_SPEC( "LinearSolverIterations", mSolNumIters );

    // compute exact residuals
    Matrix< DDRMat > tRelativeResidualNorm = mLinearSystem->compute_residual_of_linear_system();

    for ( uint i = 0; i < tRelativeResidualNorm.numel(); i++ )
    {
        MORIS_LOG_SPEC( "LinearResidualNorm", tRelativeResidualNorm( i ) );
    }

    // Delete solver
    delete mAztecSolver;
    mAztecSolver = nullptr;

    return error;
}

//----------------------------------------------------------------------------------------

void
Linear_Solver_Aztec::set_solver_internal_parameters()
{

    // Generic iterative solver parameters

    // Solver Type
    if ( mParameterList.get< moris::sint >( "AZ_solver" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_solver, mParameterList.get< moris::sint >( "AZ_solver" ) );
    }

    // Set AZ_overlap
    // Determines the submatrices factored with the domain decomposition algorithms
    // Option to specify with how many rows from other processors each processor’s local submatrix is augmented.
    if ( mParameterList.get< moris::sint >( "AZ_overlap" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_overlap, mParameterList.get< moris::sint >( "AZ_overlap" ) );
    }

    // Set AZ_type_overlap
    // AZ_standard = The resulting value of an unknown is determined by the processor owning that unknown. Information from other processors about that unknown is discarded.
    if ( mParameterList.get< moris::sint >( "AZ_type_overlap" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_type_overlap, mParameterList.get< moris::sint >( "AZ_type_overlap" ) );
    }

    // Set AZ_reorder
    // Determines whether RCM reordering will be done in conjunction with domain decomposition incomplete factorizations.
    // Option to enable (=1) or disable (=0) the Reverse Cuthill–McKee (RCM) algorithm to reorder system equations for smaller bandwidth
    if ( mParameterList.get< moris::sint >( "AZ_reorder" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_reorder, mParameterList.get< moris::sint >( "AZ_reorder" ) );
    }

    // Set AZ_aux_vec
    // AZ_resid = r_tilde is set to the initial residual vector
    // mAztecSolver.SetAztecOption ( AZ_aux_vec, AZ_resid );

    // GMRES specific solver parameters
    // Set AZ_kspace
    // Krylov subspace size for restarted GMRES
    // Setting mKrylovSpace larger improves the robustness, decreases iteration count, but increases memory consumption. For very difficult problems, set it equal to the maximum number of iterations.
    if ( mParameterList.get< moris::sint >( "AZ_kspace" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_kspace, mParameterList.get< moris::sint >( "AZ_kspace" ) );
    }

    // Set AZ_orthog
    if ( mParameterList.get< moris::sint >( "AZ_orthog" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_orthog, mParameterList.get< moris::sint >( "AZ_orthog" ) );
    }

    // Set AZ_rthresh
    // Parameter used to modify the relative magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    if ( mParameterList.get< moris::real >( "AZ_rthresh" ) != -1.0 )
    {
        mAztecSolver->SetAztecParam( AZ_rthresh, mParameterList.get< moris::real >( "AZ_rthresh" ) );
    }

    // Set AZ_athresh
    // Parameter used to modify the absolute magnitude of the diagonal entries of the matrix that is used to compute any of the incomplete factorization preconditioners
    if ( mParameterList.get< moris::real >( "AZ_athresh" ) != -1.0 )
    {
        mAztecSolver->SetAztecParam( AZ_athresh, mParameterList.get< moris::real >( "AZ_athresh" ) );
    }

    //---------------------------------------------------------------------------------------------------------------
    // Set AZ_conv criteria
    if ( mParameterList.get< moris::sint >( "AZ_conv" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_conv, mParameterList.get< moris::sint >( "AZ_conv" ) );
    }

    // Set AZ_diagnostics
    if ( mParameterList.get< moris::sint >( "AZ_diagnostics" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_diagnostics, mParameterList.get< moris::sint >( "AZ_diagnostics" ) );
    }

    // Set AZ_output
    if ( mParameterList.get< moris::sint >( "AZ_output" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_output, mParameterList.get< int >( "AZ_output" ) );
    }

    // Set if preconditioner is recalculated
    if ( mParameterList.get< moris::sint >( "AZ_pre_calc" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_pre_calc, mParameterList.get< moris::sint >( "AZ_pre_calc" ) );
    }

    // Set if preconditioner is recalculated
    if ( mParameterList.get< moris::sint >( "AZ_keep_info" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_keep_info, mParameterList.get< moris::sint >( "AZ_keep_info" ) );
    }

    // Determine which preconditioner is used
    if ( mParameterList.get< moris::sint >( "AZ_precond" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_precond, mParameterList.get< moris::sint >( "AZ_precond" ) );
    }

    // Set preconditioner subdomain solve - direct solve or incomplete
    if ( mParameterList.get< moris::sint >( "AZ_subdomain_solve" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_subdomain_solve, mParameterList.get< moris::sint >( "AZ_subdomain_solve" ) );
    }

    // Set preconditioner polynomial order
    if ( mParameterList.get< moris::sint >( "AZ_poly_ord" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_poly_ord, mParameterList.get< moris::sint >( "AZ_poly_ord" ) );
    }

    // Set drop tolerance - for LU, ILUT
    if ( mParameterList.get< moris::real >( "AZ_drop" ) != -1.0 )
    {
        mAztecSolver->SetAztecParam( AZ_drop, mParameterList.get< moris::real >( "AZ_drop" ) );
    }

    // Set level of graph fill in - for ilu(k), icc(k), bilu(k)
    if ( mParameterList.get< moris::sint >( "AZ_graph_fill" ) != INT_MAX )
    {
        mAztecSolver->SetAztecOption( AZ_graph_fill, mParameterList.get< moris::sint >( "AZ_graph_fill" ) );
    }

    // Set Damping or relaxation parameter used for RILU
    if ( mParameterList.get< moris::real >( "AZ_omega" ) != -1.0 )
    {
        mAztecSolver->SetAztecParam( AZ_omega, mParameterList.get< moris::real >( "AZ_omega" ) );
    }

    // Set ilut fill
    if ( mParameterList.get< moris::real >( "AZ_ilut_fill" ) != -1.0 )
    {
        mAztecSolver->SetAztecParam( AZ_ilut_fill, mParameterList.get< moris::real >( "AZ_ilut_fill" ) );
    }
}
