/*
 * cl_DLA_Linear_Solver_Amesos.cpp
 *
 *  Created on: May 16, 2018
 *      Author: schmidt
 */
#include "cl_DLA_Linear_Solver_Amesos.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_Vector.hpp"

#include "cl_Tracer.hpp"

using namespace moris;
using namespace dla;

Linear_Solver_Amesos::Linear_Solver_Amesos( Linear_Problem * aLinearSystem ) //: mAmesosSolver( NULL ),
                                                                         //  mAmesosFactory()
{
    // boolean for symbolic factorization after first solve
    mIsPastFirstSolve = false;

    // Set chosen solver options
    this->set_solver_parameters();
}

Linear_Solver_Amesos::~Linear_Solver_Amesos()
{
//    delete mAmesosSolver;
}

void Linear_Solver_Amesos::set_solver_parameters()
{
    // ASSIGN DEFAULT PARAMETER VALUES
    // Amesos 2.0 Reference Guide, SANDIA REPORT, SAND2004-4820, https://trilinos.org/oldsite/packages/amesos/AmesosReferenceGuide.pdf

    // Set Amesos solver type
    // options are "Amesos_Klu", "Amesos_Superlu", "Amesos_Umfpack", "Amesos_Dscpack",
    // "Amesos_Superludist", "Amesos_Mumps", "Amesos_Scalapack", "Amesos_Pardiso"
//    mParameterList.insert( "solver_type" ,  "Amesos_Pardiso" );

    // set AZ_output options
    // options are true, false
//    mParameterList.insert( "output" , false );

    // set symbolic factorization
    // options are true, false
//    mParameterList.insert( "symbolic_factorization" , false );
}

moris::sint Linear_Solver_Amesos::solve_linear_system()
{
    sint error = 0;
    moris::real startSolTime     = 0.0;
    moris::real startSymFactTime = 0.0;
    moris::real startNumFactTime = 0.0;

    // Set all Aztec options
    this->set_solver_internal_parameters();

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
    const moris::real endNumFactTime = ( mAmesosSolver->NumNumericFact()  > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total numeric factorization time"  ) : 0.0;

    mSolTime     = endSolTime     - startSolTime;
    mSymFactTime = endSymFactTime - startSymFactTime;
    mNumFactTime = endNumFactTime - startNumFactTime;

    return error;
}

moris::sint Linear_Solver_Amesos::solve_linear_system( Linear_Problem * aLinearSystem, const moris::sint aIter )
{
    Tracer tTracer(EntityBase::LinearSolver, EntityType::Amesos, EntityAction::Solve);

    mLinearSystem = aLinearSystem;
    mAmesosSolver = mAmesosFactory.Create( "Amesos_Pardiso", *mLinearSystem->get_linear_system_epetra() );

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
        error = mAmesosSolver->SymbolicFactorization();

        MORIS_ERROR( error == 0, "SYMBOLIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve", error );
//    }

    // Perform numeric factorization
    error = mAmesosSolver->NumericFactorization();
    MORIS_ERROR( error == 0, "NUMERIC FACTORIZATION in Linear Solver Trilinos Amesos returned an error %i. Exiting linear solve", error );

    // Solve linear system
    error = mAmesosSolver->Solve();
    MORIS_ERROR( error == 0, "Error in solving linear system with Amesos" );

    mIsPastFirstSolve = true;

    // Get chosen times from solver
    mAmesosSolver->GetTiming( timingsList );

    const moris::real endSolTime     = Teuchos::getParameter< moris::real >( timingsList, "Total solve time" );
    const moris::real endSymFactTime = ( mAmesosSolver->NumSymbolicFact() > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total symbolic factorization time" ) : 0.0;
    const moris::real endNumFactTime = ( mAmesosSolver->NumNumericFact()  > 0 ) ? Teuchos::getParameter< moris::real >( timingsList, "Total numeric factorization time"  ) : 0.0;

    mSolTime     = endSolTime     - startSolTime;
    mSymFactTime = endSymFactTime - startSymFactTime;
    mNumFactTime = endNumFactTime - startNumFactTime;

    return error;
}

void Linear_Solver_Amesos::set_solver_internal_parameters()
{
    // Set Amesos solver type
//    mAmesosSolver = mAmesosFactory.Create( mParameterList.get< const char * >( "solver_type" ), mEpetraProblem );
//    mAmesosSolver = mAmesosFactory.Create( "Amesos_Klu", *mLinearSystem->get_linear_system_epetra() );

    // Initialize parameter list
//    Teuchos::ParameterList params;

    // Set output options
//    params.set("PrintTiming"        , mParameterList.get< bool >( "output" ));
//    params.set("PrintStatus"        , mParameterList.get< bool >( "output" ));
//    params.set("ComputeVectorNorms" , mParameterList.get< bool >( "output" ));
//    params.set("ComputeTrueResidual", mParameterList.get< bool >( "output" ));

    // Allows non contiguous indexing of stiffness matrix dofs (XFEM Reduced system) KLU only
//    params.set("Reindex",(bool) (true));

    // Create a parameter sublist for PARDISO solver
//    const char* tSolverType1 = "Amesos_Pardiso";
////    if ( mParameterList.get< const char * >( "solver_type" ) == tSolverType1 )
//    if ( "Amesos_Pardiso" == tSolverType1 )
    {
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
    }

//    // Create a parameter sub list for MUMPS solver
//    const char* tSolverType2 = "Amesos_Mumps";
////    if (mParameterList.get< const char * >( "solver_type" ) == tSolverType2 )
//    if ( "Amesos_Mumps" == tSolverType2 )
//    {
//        // Increase memory allocation to 200% at a time (default = 20%)
//        params.sublist( "mumps" ).set( "ICNTL(14)",static_cast<int>(200) );
//    }

    // Set AMESOS solver parameter list
//    mAmesosSolver->SetParameters( params );
}
