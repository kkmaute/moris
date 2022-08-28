/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Factory.cpp
 *
 */

#include "cl_TSA_Time_Solver_Factory.hpp"

#include "cl_TSA_Monolithic_Time_Solver.hpp"
#include "cl_TSA_Staggered_Time_Solver.hpp"

using namespace moris;
using namespace tsa;

Time_Solver_Factory::Time_Solver_Factory()
{}

Time_Solver_Factory::~Time_Solver_Factory()
{}

std::shared_ptr< Time_Solver_Algorithm > Time_Solver_Factory::create_time_solver( const enum TimeSolverType aTimeSolverType )
{
    std::shared_ptr< Time_Solver_Algorithm > tTimeSolver = nullptr;

    switch( aTimeSolverType )
    {
    case ( TimeSolverType::MONOLITHIC ):
        tTimeSolver = std::make_shared< Monolithic_Time_Solver >();
        break;
    case ( TimeSolverType::STAGGERED ):
        tTimeSolver = std::make_shared<  Staggered_Time_Solver >();
        break;
    default:
        MORIS_ERROR( false, "No solver type specified" );
        break;
    }

    return tTimeSolver;
}

std::shared_ptr< Time_Solver_Algorithm > Time_Solver_Factory::create_time_solver( const enum TimeSolverType aTimeSolverType,
                                                                                  const ParameterList       aParameterlist )
{
    std::shared_ptr< Time_Solver_Algorithm > tTimeSolver = nullptr;

    switch( aTimeSolverType )
    {
    case ( TimeSolverType::MONOLITHIC ):
        tTimeSolver = std::make_shared< Monolithic_Time_Solver >( aParameterlist );
        break;
    case ( TimeSolverType::STAGGERED ):
        tTimeSolver = std::make_shared<  Staggered_Time_Solver >( aParameterlist );
        break;
    default:
        MORIS_ERROR( false, "No solver type specified" );
        break;
    }

    return tTimeSolver;
}

