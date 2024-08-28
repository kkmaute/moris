/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Solver_Factory.cpp
 *
 */

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Linear_Solver_Amesos.hpp"
#include "cl_DLA_Linear_Solver_Belos.hpp"
#include "cl_DLA_Linear_Solver_ML.hpp"

#include "cl_DLA_Linear_System_Trilinos.hpp"

#ifdef MORIS_HAVE_PETSC
#include "cl_DLA_Linear_System_PETSc.hpp"
#include "cl_DLA_Linear_Solver_PETSc.hpp"
#ifdef MORIS_HAVE_SLEPC
#include "cl_DLA_Eigen_Solver_SLEPc.hpp"
#endif  
#endif

#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"

using namespace moris;
using namespace dla;

//------------------------------------------------------------------------------

Solver_Factory::Solver_Factory()
{
}

//------------------------------------------------------------------------------

Solver_Factory::~Solver_Factory()
{
}

//------------------------------------------------------------------------------
Preconditioner*
Solver_Factory::create_preconditioner( const Parameter_List& aParameterList )
{
    switch ( aParameterList.get< sol::PreconditionerType >( "Preconditioner_Implementation" ) )
    {
        case ( sol::PreconditionerType::NONE ):
            return nullptr;
        case ( sol::PreconditionerType::IFPACK ):
        case ( sol::PreconditionerType::ML ):
            return new Preconditioner_Trilinos( aParameterList );
        case ( sol::PreconditionerType::PETSC ):
            return new Preconditioner_PETSc( aParameterList );
        default:
            MORIS_ERROR( false, "No solver type specified" );
            return nullptr;
    }
}

//------------------------------------------------------------------------------

std::shared_ptr< Linear_Solver_Algorithm >
Solver_Factory::create_solver(
        const Parameter_List& aParameterlist )
{
    std::shared_ptr< Linear_Solver_Algorithm > tLinSol;

    switch ( aParameterlist.get< sol::SolverType >( "Solver_Implementation" ) )
    {
        case ( sol::SolverType::AZTEC_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Aztec >( aParameterlist );
            break;
        case ( sol::SolverType::AMESOS_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Amesos >( aParameterlist );
            break;
        case ( sol::SolverType::BELOS_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Belos >( aParameterlist );
            break;
            //    case ( sol::SolverType::AMESOS2_IMPL ):
            //        tLinSol = std::make_shared< Linear_Solver_Amesos2 >( aParameterlist );
            //        break;
        case ( sol::SolverType::PETSC ):
#ifdef MORIS_HAVE_PETSC
            tLinSol = std::make_shared< Linear_Solver_PETSc >( aParameterlist );
#else
            MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
            break;
        case ( sol::SolverType::EIGEN_SOLVER ):
            tLinSol = std::make_shared< Eigen_Solver >( &aParameterlist );
            break;
        case ( sol::SolverType::ML ):
            tLinSol = std::make_shared< Linear_Solver_ML >( aParameterlist );
            break;
        case ( sol::SolverType::SLEPC_SOLVER ):
#ifdef MORIS_HAVE_SLEPC
            tLinSol = std::make_shared< Eigen_Solver_SLEPc >(aParameterlist);
#else
            MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
            break;
        default:
            MORIS_ERROR( false, "No solver type specified" );
            break;
    }
    return tLinSol;
}

//------------------------------------------------------------------------------

Linear_Problem*
Solver_Factory::create_linear_system(
        moris::Solver_Interface* aSolverInterface,
        sol::SOL_Warehouse*      aSolverWarehouse,
        sol::Dist_Map*           aMap,
        sol::Dist_Map*           aFullMap,
        const enum sol::MapType  aLinSysType,
        const bool               aNotCreatedByNonLinSolver )
{
    Linear_Problem* tLinSys;

    switch ( aLinSysType )
    {
        case ( sol::MapType::Epetra ):
            tLinSys = new Linear_System_Trilinos( aSolverInterface, aSolverWarehouse, aMap, aFullMap );
            break;
        case ( sol::MapType::Petsc ):
#ifdef MORIS_HAVE_PETSC
            tLinSys = new Linear_System_PETSc(
                    aSolverInterface,
                    aSolverWarehouse,
                    aMap,
                    aFullMap,
                    aNotCreatedByNonLinSolver );
#else
            MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
            break;
        default:
            MORIS_ERROR( false, "Solver_Factory::create_linear_system: No solver type specified" );
            break;
    }

    return tLinSys;
}

//------------------------------------------------------------------------------

Linear_Problem*
Solver_Factory::create_linear_system(
        moris::Solver_Interface* aSolverInterface,
        const enum sol::MapType  aLinSysType,
        const bool               aNotCreatedByNonLinSolver )
{
    Linear_Problem* tLinSys = nullptr;

    switch ( aLinSysType )
    {
        case ( sol::MapType::Epetra ):
            tLinSys = new Linear_System_Trilinos( aSolverInterface );
            break;
        case ( sol::MapType::Petsc ):
#ifdef MORIS_HAVE_PETSC
            tLinSys = new Linear_System_PETSc( aSolverInterface, aNotCreatedByNonLinSolver );
#else
            MORIS_ERROR( false, "MORIS is configured with out PETSC support." );
#endif
            break;
        default:
            MORIS_ERROR( false, "Solver_Factory::create_linear_system: No solver type specified" );
            break;
    }

    return tLinSys;
}
