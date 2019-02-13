/*
 * cl_NLA_Nonlinear_Solver_Factory.cpp
 *
 *  Created on: Apr 10, 2018
 *      Author: schmidt
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

Time_Solver * Time_Solver_Factory::create_time_solver( const enum TimeSolverType aTimeSolverType )
{
     Time_Solver * tTimeSolver = nullptr;

    switch( aTimeSolverType )
    {
    case ( TimeSolverType::MONOLITHIC ):
        tTimeSolver = new Monolithic_Time_Solver();
        break;
    case ( TimeSolverType::STAGGERED ):
        tTimeSolver = new Staggered_Time_Solver ();
        break;
    default:
        MORIS_ERROR( false, "No solver type specified" );
        break;
    }

    return tTimeSolver;
}

