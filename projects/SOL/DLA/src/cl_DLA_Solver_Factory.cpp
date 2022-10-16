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
#include "cl_DLA_Linear_Solver_Amesos2.hpp"
#include "cl_DLA_Linear_Solver_Belos.hpp"
#include "cl_DLA_Linear_Solver_PETSc.hpp"

#include "cl_DLA_Linear_System_Trilinos.hpp"
#include "cl_DLA_Linear_System_PETSc.hpp"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

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

std::shared_ptr< Linear_Solver_Algorithm >
Solver_Factory::create_solver( const enum sol::SolverType aSolverType )
{
    std::shared_ptr< Linear_Solver_Algorithm > tLinSol;

    switch ( aSolverType )
    {
        case ( sol::SolverType::AZTEC_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Aztec >();
            break;
        case ( sol::SolverType::AMESOS_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Amesos >();
            break;
        case ( sol::SolverType::BELOS_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Belos >();
            break;
        case ( sol::SolverType::AMESOS2_IMPL ):
            tLinSol = std::make_shared< Linear_Solver_Amesos2 >();
            break;
        case ( sol::SolverType::PETSC ):
            tLinSol = std::make_shared< Linear_Solver_PETSc >();
            break;
        default:
            MORIS_ERROR( false, "No solver type specified" );
            break;
    }
    return tLinSol;
}

//------------------------------------------------------------------------------

std::shared_ptr< Linear_Solver_Algorithm >
Solver_Factory::create_solver(
        const enum sol::SolverType aSolverType,
        const ParameterList        aParameterlist )
{
    std::shared_ptr< Linear_Solver_Algorithm > tLinSol;

    switch ( aSolverType )
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
            tLinSol = std::make_shared< Linear_Solver_PETSc >( aParameterlist );
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
            tLinSys = new Linear_System_PETSc(
                    aSolverInterface,
                    aSolverWarehouse,
                    aMap,
                    aFullMap,
                    aNotCreatedByNonLinSolver );
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
            tLinSys = new Linear_System_PETSc( aSolverInterface, aNotCreatedByNonLinSolver );
            break;
        default:
            MORIS_ERROR( false, "Solver_Factory::create_linear_system: No solver type specified" );
            break;
    }

    return tLinSys;
}
