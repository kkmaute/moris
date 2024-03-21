/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Belos.cpp
 *
 */

#include "cl_DLA_Linear_Solver_Belos.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

#include "cl_DLA_Preconditioner_Trilinos.hpp"

#include "cl_Tracer.hpp"
#include "cl_Logger.hpp"

#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>

#include <BelosSolverManager.hpp>
#include "BelosSolverFactory.hpp"

#include "Teuchos_ParameterList.hpp"

#include "BelosEpetraAdapter.hpp"

#include "Ifpack.h"
#include "Ifpack_AdditiveSchwarz.h"

using namespace moris;
using namespace dla;

Linear_Solver_Belos::Linear_Solver_Belos( const moris::Parameter_List& aParameterlist )
        : Linear_Solver_Algorithm_Trilinos( aParameterlist )
{
}

//---------------------------------------------------------------------------------------------------

Linear_Solver_Belos::Linear_Solver_Belos( Linear_Problem* aLinearSystem )
        : Linear_Solver_Algorithm_Trilinos( prm::create_linear_algorithm_parameter_list_belos() )
{
}

//---------------------------------------------------------------------------------------------------

Linear_Solver_Belos::~Linear_Solver_Belos()
{
}

//---------------------------------------------------------------------------------------------------

moris::sint
Linear_Solver_Belos::solve_linear_system(
        Linear_Problem*   aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "Belos", "Solve" );

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ParameterList;
    using Teuchos::parameterList;
    using Belos::SolverFactory;

    // set Belos solver options
    this->set_solver_internal_parameters();

    // set linear system
    mLinearSystem = aLinearSystem;

    // mPreconditioner->initialize( &mParameterList, mLinearSystem );
    mPreconditioner->build( aLinearSystem, aIter );

    MORIS_ERROR( mPreconditioner->exists(),
            "Linear_Solver_Belos::solve_linear_system - No preconditioner has been defined.\n" );

    RCP< Belos::EpetraPrecOp > belosPrec =
            rcp( new Belos::EpetraPrecOp( mPreconditioner->get_operator() ) );

    // get operator, solution and Rhs vectors
    RCP< Epetra_CrsMatrix > A =
            rcp( dynamic_cast< Epetra_CrsMatrix* >( aLinearSystem->get_matrix()->get_matrix() ), false );
    RCP< Epetra_MultiVector > X =
            rcp( dynamic_cast< Vector_Epetra* >( aLinearSystem->get_free_solver_LHS() )->get_epetra_vector(), false );
    RCP< Epetra_MultiVector > B =
            rcp( dynamic_cast< Vector_Epetra* >( aLinearSystem->get_solver_RHS() )->get_epetra_vector(), false );

    // create linear problem
    RCP< Belos::LinearProblem< double, Epetra_MultiVector, Epetra_Operator > > problem =
            rcp( new Belos::LinearProblem< double, Epetra_MultiVector, Epetra_Operator >( A, X, B ) );

    // set either left or right preconditioner
    if ( mParameterList.get< std::string >( "Left-right Preconditioner" ) == "left" )
    {
        problem->setLeftPrec( belosPrec );
    }
    else
    {
        problem->setRightPrec( belosPrec );
    }

    // check problem set
    MORIS_ERROR( problem->setProblem(),
            "Linear_Solver_Belos::solve_linear_system - LinearProblem is not correctly set up.\n" );

    // Create iterative solver.
    SolverFactory< double, Epetra_MultiVector, Epetra_Operator > factory;

    RCP< Belos::SolverManager< double, Epetra_MultiVector, Epetra_Operator > > solver =
            factory.create(
                    mParameterList.get< std::string >( "Solver Type" ),
                    mMyPl );

    // Tell the solver what problem you want to solve.
    solver->setProblem( problem );

    // Solve problem
    Belos::ReturnType tSolverConvergence = solver->solve();

    MORIS_LOG_SPEC( "IterativeSolverConverged", tSolverConvergence );

    // Ask the solver how many iterations the last solve() took.
    MORIS_LOG_SPEC( "LinearSolverIterations", solver->getNumIters() );

    // Get solution tolerance across all RHS
    MORIS_LOG_SPEC( "LinearResidualNorm_All_RHS", solver->achievedTol() );

    // compute exact residuals
    Matrix< DDRMat > tRelativeResidualNorm = aLinearSystem->compute_residual_of_linear_system();

    for ( uint i = 0; i < tRelativeResidualNorm.numel(); i++ )
    {
        MORIS_LOG_SPEC( "LinearResidualNorm_RHS_" + std::to_string( i ), tRelativeResidualNorm( i ) );
    }

    // return solver status
    return 0;
}

//---------------------------------------------------------------------------------------------------

void
Linear_Solver_Belos::set_solver_internal_parameters()
{
    mMyPl = Teuchos::parameterList();

    if ( mParameterList.get< moris::sint >( "Verbosity" ) != INT_MAX )
    {
        mMyPl->set( "Verbosity", mParameterList.get< moris::sint >( "Verbosity" ) );
    }

    if ( mParameterList.get< moris::sint >( "Num Blocks" ) != INT_MAX )
    {
        mMyPl->set( "Num Blocks", mParameterList.get< moris::sint >( "Num Blocks" ) );
    }

    if ( mParameterList.get< moris::sint >( "Block Size" ) != INT_MAX )
    {
        mMyPl->set( "Block Size", mParameterList.get< moris::sint >( "Block Size" ) );
    }

    if ( mParameterList.get< moris::sint >( "Maximum Iterations" ) != INT_MAX )
    {
        mMyPl->set( "Maximum Iterations", mParameterList.get< moris::sint >( "Maximum Iterations" ) );
    }

    if ( mParameterList.get< moris::real >( "Convergence Tolerance" ) != 1e-08 )
    {
        mMyPl->set( "Convergence Tolerance", mParameterList.get< moris::real >( "Convergence Tolerance" ) );
    }

    if ( mParameterList.get< moris::sint >( "Output Frequency" ) != -1 )
    {
        mMyPl->set( "Output Frequency", mParameterList.get< moris::real >( "Output Frequency" ) );
    }
}
