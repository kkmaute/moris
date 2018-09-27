/*
 * cl_Linear_Solver_Amesos2.cpp
 *
 *  Created on: May 16, 2018
 *      Author: schmidt
 */

#include <string.h>

#include "cl_Linear_Solver_Amesos2.hpp"

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

using namespace moris;

Linear_Solver_Amesos2::Linear_Solver_Amesos2( Solver_Input*   aInput ) : Linear_Solver_Trilinos ( aInput )
{
    // boolean for symbolic factorization after first solve
    mIsPastFirstSolve = false;

    using Tpetra::global_size_t;
    using Teuchos::tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;

    // Teuchos smart pointer to Matrix, RHS and LHS

    mAMatrix   = rcp( dynamic_cast< Epetra_CrsMatrix* > ( mEpetraProblem.GetMatrix()) );
    mRHS       = rcp( mEpetraProblem.GetRHS() );
    mLHS       = rcp( mEpetraProblem.GetLHS() );

    // Set chosen solver options
    this->set_solver_parameters();
}

Linear_Solver_Amesos2::~Linear_Solver_Amesos2()
{
    //    Delete( mAmesos2Solver);
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
    //int error = 0;

    // Set all Aztec options
    this->set_solver_internal_parameters();

    // Get timing info
    Teuchos::ParameterList timingsList;

    if ( !mParameterList.get< bool >( "symbolic_factorization" ) || !mIsPastFirstSolve )
    {
        mAmesos2Solver->symbolicFactorization();
    }

    mAmesos2Solver->numericFactorization();

    mAmesos2Solver->solve();

    mIsPastFirstSolve = true;

    mAmesos2Solver->getTiming(timingsList);

    mSolTime     = 0;

    return 0;
}

void Linear_Solver_Amesos2::set_solver_internal_parameters()
{
    Teuchos::ParameterList params( "Amesos2" );

    mAmesos2Solver = Amesos2::create<Epetra_CrsMatrix,Epetra_MultiVector>( mParameterList.get< const char * >( "solver_type" ), mAMatrix, mLHS, mRHS );

    const char* tSolverType1 = "pardiso_mkl";
    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType1 )
    {
        params.set<bool>( "Transpose", false );
        params.sublist( "PARDISOMKL" ).set< std::string >( "IPARM(12)", "0" );
    }

    const char* tSolverType2 = "mumps";
    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType2 )
    {

    }

    mAmesos2Solver->setParameters( rcpFromRef( params ) );
}

