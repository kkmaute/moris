/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Amesos2.cpp
 *
 */

#include <string.h>

#include "cl_DLA_Linear_Solver_Amesos2.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_Amesos2::Linear_Solver_Amesos2( Solver_Interface *   aInput )
{
//    // boolean for symbolic factorization after first solve
//    mIsPastFirstSolve = false;
//
//    using Tpetra::global_size_t;
//    using Teuchos::tuple;
//    using Teuchos::RCP;
//    using Teuchos::rcp;
//
//    // Teuchos smart pointer to Matrix, RHS and LHS
//
//    mAMatrix   = rcp( dynamic_cast< Epetra_CrsMatrix* > ( mEpetraProblem.GetMatrix()) );
//    mRHS       = rcp( mEpetraProblem.GetRHS() );
//    mLHS       = rcp( mEpetraProblem.GetLHS() );
//
//    // Set chosen solver options
//    this->set_solver_parameters();
}

Linear_Solver_Amesos2::Linear_Solver_Amesos2()
{

}

Linear_Solver_Amesos2::~Linear_Solver_Amesos2()
{
    // Set chosen solver options
    this->set_solver_parameters();
}

void Linear_Solver_Amesos2::set_solver_parameters()
{
    // ASSIGN DEFAULT PARAMETER VALUES
    // Amesos 2.0 Reference Guide, SANDIA REPORT, SAND2004-4820, https://trilinos.org/oldsite/packages/amesos/AmesosReferenceGuide.pdf

    // Set Amesos solver type
    // options are "mumps", "pardiso_mkl"
    mParameterList.insert( "solver_type" ,  "pardiso_mkl" );

    // set AZ_output options
    // options are true, false
    mParameterList.insert( "output" , false );

    // set symbolic factorization
    // options are true, false
    mParameterList.insert( "symbolic_factorization" , false );
}

moris::sint Linear_Solver_Amesos2::solve_linear_system()
{
//    //int error = 0;
//
//    // Set all Aztec options
//    this->set_solver_internal_parameters();
//
//    // Get timing info
//    Teuchos::ParameterList timingsList;
//
//    if ( !mParameterList.get< bool >( "symbolic_factorization" ) || !mIsPastFirstSolve )
//    {
//        mAmesos2Solver->symbolicFactorization();
//    }
//
//    mAmesos2Solver->numericFactorization();
//
//    mAmesos2Solver->solve();
//
//    mIsPastFirstSolve = true;
//
//    mAmesos2Solver->getTiming(timingsList);
//
//    mSolTime     = 0;
//
    return 0;
}

moris::sint Linear_Solver_Amesos2::solve_linear_system(       Linear_Problem * aLinearSystem,
                                                        const moris::sint     aIter )
{
//    this->set_solver_internal_parameters();
//    mLinearSystem = aLinearSystem;
//
//    //    //int error = 0;
//
//    using Teuchos::RCP;
//    using Teuchos::rcp;
//    using Teuchos::ParameterList;
//    using Teuchos::parameterList;
//    using Belos::SolverFactory;
//
//
////    RCP<Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> > problem
////                   = rcp (new Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> (rcp( dynamic_cast< Epetra_CrsMatrix* > ( aLinearSystem->get_matrix()->get_matrix() ), false ),
////                                                                                               rcp( aLinearSystem->get_free_solver_LHS()->get_vector(), false ),
////                                                                                               rcp( aLinearSystem->get_solver_RHS()->get_vector(), false ) ) );
//
//    tAMatrix   = Teuchos::rcp( dynamic_cast< Epetra_CrsMatrix* > ( aLinearSystem->get_matrix()->get_matrix() ), false );
//    tRHS       = Teuchos::rcp( aLinearSystem->get_free_solver_LHS()->get_vector()                             , false );
//    tLHS       = Teuchos::rcp( aLinearSystem->get_solver_RHS()     ->get_vector()                             , false );
//
//    mAmesos2Solver = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>( "pardiso_mkl", tAMatrix, tLHS, tRHS );
//
//    // Get timing info
//    Teuchos::ParameterList timingsList;
//
//    if ( !mParameterList.get< bool >( "symbolic_factorization" ) || !mIsPastFirstSolve )
//    {
//        mAmesos2Solver->symbolicFactorization();
//    }
//
//    mAmesos2Solver->numericFactorization();
//
//    mAmesos2Solver->solve();
//
//    mIsPastFirstSolve = true;
//
//    mAmesos2Solver->getTiming(timingsList);
//
//    mSolTime     = 0;
//
    return 0;
}

void Linear_Solver_Amesos2::set_solver_internal_parameters()
{
//    Teuchos::ParameterList params( "Amesos2" );
//
//    mAmesos2Solver = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>( mParameterList.get< const char * >( "solver_type" ), mAMatrix, mLHS, mRHS );
//
//    const char* tSolverType1 = "pardiso_mkl";
//    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType1 )
//    {
//        params.set<bool>( "Transpose", false );
//        params.sublist( "PARDISOMKL" ).set< std::string >( "IPARM(12)", "0" );
//    }
//
//    const char* tSolverType2 = "mumps";
//    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType2 )
//    {
//
//    }
//
//    mAmesos2Solver->setParameters( rcpFromRef( params ) );
}

