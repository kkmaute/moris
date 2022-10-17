/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Amesos.cpp
 *
 */

#include "cl_DLA_Linear_Solver_Amesos.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_SOL_Dist_Vector.hpp"
#include "cl_SOL_Dist_Matrix.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

#include "Amesos_Umfpack.h"

#include "cl_Tracer.hpp"

using namespace moris;
using namespace dla;

//-----------------------------------------------------------------------------

Linear_Solver_Amesos::Linear_Solver_Amesos(
        const moris::ParameterList aParameterlist )
        : Linear_Solver_Algorithm( aParameterlist )
{
    // boolean for symbolic factorization after first solve
    mIsPastFirstSolve = false;
}

//-----------------------------------------------------------------------------

Linear_Solver_Amesos::Linear_Solver_Amesos( Linear_Problem* aLinearSystem )
{
    // boolean for symbolic factorization after first solve
    mIsPastFirstSolve = false;

    // Set chosen solver options
    this->set_solver_parameters();
}

Linear_Solver_Amesos::Linear_Solver_Amesos()
{
    // Set chosen solver options
    this->set_solver_parameters();
}

//-----------------------------------------------------------------------------

Linear_Solver_Amesos::~Linear_Solver_Amesos()
{
    delete mAmesosSolver;
    mAmesosSolver = nullptr;
}

//-----------------------------------------------------------------------------

void
Linear_Solver_Amesos::set_solver_parameters()
{
    mParameterList = prm::create_linear_algorithm_parameter_list_amesos();
}

//-----------------------------------------------------------------------------

moris::sint
Linear_Solver_Amesos::solve_linear_system()
{
    Tracer tTracer( "LinearAlgorithm", "Amesos", "Solve" );

    sint error = 0;

    moris::real startSolTime     = 0.0;
    moris::real startSymFactTime = 0.0;
    moris::real startNumFactTime = 0.0;

    // Get timing info
    Teuchos::ParameterList timingsList;

    // Perform symbolic factorization
    if ( !mParameterList.get< bool >( "symbolic_factorization" ) || !mIsPastFirstSolve )
    {
        error = mAmesosSolver->SymbolicFactorization();
        MORIS_ERROR( error == 0, "SYMBOLIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve" );
    }

    // Perform numeric factorization
    error = mAmesosSolver->NumericFactorization();
    MORIS_ERROR( error == 0, "NUMERIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve" );

    // Solve linear system
    error = mAmesosSolver->Solve();
    MORIS_ERROR( error == 0, "Error in solving linear system with Amesos" );

    mIsPastFirstSolve = true;

    // Get chosen times from solver
    mAmesosSolver->GetTiming( timingsList );

    const moris::real endSolTime     = Teuchos::getParameter< moris::real >( timingsList, "Total solve time" );
    const moris::real endSymFactTime = ( mAmesosSolver->NumSymbolicFact() > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total symbolic factorization time" ) : 0.0;
    const moris::real endNumFactTime = ( mAmesosSolver->NumNumericFact() > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total numeric factorization time" ) : 0.0;

    mSolTime     = endSolTime - startSolTime;
    mSymFactTime = endSymFactTime - startSymFactTime;
    mNumFactTime = endNumFactTime - startNumFactTime;

    delete mAmesosSolver;
    mAmesosSolver = nullptr;

    return error;
}

//-----------------------------------------------------------------------------

moris::sint
Linear_Solver_Amesos::solve_linear_system(
        Linear_Problem*   aLinearSystem,
        const moris::sint aIter )
{
    Tracer tTracer( "LinearSolver", "Amesos", "Solve" );

    mLinearSystem = aLinearSystem;

    mEpetraProblem.SetOperator( aLinearSystem->get_matrix()->get_matrix() );
    mEpetraProblem.SetRHS( dynamic_cast< Vector_Epetra* >( aLinearSystem->get_solver_RHS() )->get_epetra_vector() );
    mEpetraProblem.SetLHS( dynamic_cast< Vector_Epetra* >( aLinearSystem->get_free_solver_LHS() )->get_epetra_vector() );

    Amesos tAmesosFactory;

    std::string tSolverType = mParameterList.get< std::string >( "Solver_Type" );

    mAmesosSolver = tAmesosFactory.Create( tSolverType, mEpetraProblem );

    MORIS_ERROR( mAmesosSolver,
            "Linear_Solver_Amesos::solve_linear_system - solver not implemented" );

    sint error = 0;

    moris::real startSolTime     = 0.0;
    moris::real startSymFactTime = 0.0;
    moris::real startNumFactTime = 0.0;

    // Set all Aztec options
    this->set_solver_internal_parameters();

    // Get timing info
    Teuchos::ParameterList timingsList;

    // Perform symbolic factorization
    //    if ( !mParameterList.get< bool >( "symbolic_factorization" ) || !mIsPastFirstSolve )
    //    {
    // FIXME Symbilic factorization only if matrix graph changes
    error = mAmesosSolver->SymbolicFactorization();

    MORIS_ERROR( error == 0, "SYMBOLIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve", error );
    //    }

    // Perform numeric factorization // FIXME only if values of matrix change
    error = mAmesosSolver->NumericFactorization();
    MORIS_ERROR( error == 0, "NUMERIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve", error );

    // Solve linear system
    error = mAmesosSolver->Solve();
    MORIS_ERROR( error == 0, "Error in solving linear system with Amesos" );

    mIsPastFirstSolve = true;

    // Get chosen times from solver
    mAmesosSolver->GetTiming( timingsList );

    // FIXME: should be only on demand - not by default
    if ( tSolverType == "Amesos_Umfpack" )
    {
        real tConditionNumber = reinterpret_cast< Amesos_Umfpack* >( mAmesosSolver )->GetRcond();

        MORIS_LOG_SPEC( "Condition Number of Operator: ", tConditionNumber );
    }

    const moris::real endSolTime = Teuchos::getParameter< moris::real >( timingsList, "Total solve time" );

    const moris::real endSymFactTime = ( mAmesosSolver->NumSymbolicFact() > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total symbolic factorization time" ) : 0.0;
    const moris::real endNumFactTime = ( mAmesosSolver->NumNumericFact() > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total numeric factorization time" ) : 0.0;

    mSolTime     = endSolTime - startSolTime;
    mSymFactTime = endSymFactTime - startSymFactTime;
    mNumFactTime = endNumFactTime - startNumFactTime;

    delete mAmesosSolver;
    mAmesosSolver = nullptr;

    return error;
}

//-----------------------------------------------------------------------------

void
Linear_Solver_Amesos::set_solver_internal_parameters()
{
    // Initialize parameter list
    Teuchos::ParameterList params;

    if ( mParameterList.get< moris::real >( "RcondThreshold" ) != -1.0 )
    {
        params.set( "RcondThreshold", mParameterList.get< moris::real >( "RcondThreshold" ) );
    }

    if ( mParameterList.get< moris::sint >( "OutputLevel" ) != INT_MAX )
    {
        params.set( "OutputLevel", mParameterList.get< moris::sint >( "OutputLevel" ) );
    }

    if ( mParameterList.get< moris::sint >( "DebugLevel" ) != INT_MAX )
    {
        params.set( "DebugLevel", mParameterList.get< moris::sint >( "DebugLevel" ) );
    }

    params.set( "PrintStatus", mParameterList.get< bool >( "PrintStatus" ) );
    params.set( "PrintTiming", mParameterList.get< bool >( "PrintTiming" ) );
    params.set( "ComputeVectorNorms", mParameterList.get< bool >( "ComputeVectorNorms" ) );
    params.set( "ComputeTrueResidual", mParameterList.get< bool >( "ComputeTrueResidual" ) );
    params.set( "Reindex", mParameterList.get< bool >( "Reindex" ) );
    params.set( "Refactorize", mParameterList.get< bool >( "Refactorize" ) );
    params.set( "AddZeroToDiag", mParameterList.get< bool >( "AddZeroToDiag" ) );

    std::string tSolverType = mParameterList.get< std::string >( "Solver_Type" );

    // Create a parameter sub list for MUMPS solver
    if ( tSolverType == "Amesos_Mumps" )
    {
        // Increase memory allocation to 200% at a time (default = 20%)
        params.sublist( "mumps" ).set( "ICNTL(14)", static_cast< int >( 200 ) );
    }

    // Create a parameter sublist for PARDISO solver
    //    const char* tSolverType1 = "Amesos_Pardiso";
    ////    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType1 )
    //    if ( "Amesos_Pardiso" == tSolverType1 )
    // {
    //        params.sublist( "Pardiso" ).set( "Use (non-)symmetric scaling vectors", static_cast<int>(1) );
    //        params.sublist( "Pardiso" ).set( "Use (non-)symmetric matchings"      , static_cast<int>(1) );
    //        if (aUseTranspose)
    //        {
    //            params.sublist("Pardiso").set("Solve transposed",static_cast<int>(1) );
    //        }
    //        else
    //        {
    //            params.sublist( "Pardiso" ).set( "Solve transposed", static_cast<int>(0) );
    //        }
    //}

    // Set AMESOS solver parameter list
    mAmesosSolver->SetParameters( params );
}
