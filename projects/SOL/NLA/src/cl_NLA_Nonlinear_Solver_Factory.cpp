/*
 * cl_NLA_Nonlinear_Solver_Factory.cpp
 *
 *  Created on: Apr 10, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Newton_Solver.hpp"

using namespace moris;
using namespace NLA;

Nonlinear_Solver_Factory::Nonlinear_Solver_Factory()
{}

Nonlinear_Solver_Factory::~Nonlinear_Solver_Factory()
{}

std::shared_ptr< Nonlinear_Solver > Nonlinear_Solver_Factory::create_nonlinear_solver(       Solver_Interface         * aSolverInput,
                                                                                       const enum NonlinearSolverType   aNonLinSolverType )
{
    std::shared_ptr< Nonlinear_Solver > tNonLinSys;

    switch( aNonLinSolverType )
    {
    case ( NonlinearSolverType::NEWTON_SOLVER ):
        tNonLinSys = std::make_shared< Newton_Solver >( aSolverInput );
        break;
    default:
        MORIS_ASSERT( false, "No solver type specified" );
        break;
    }

    return tNonLinSys;
}

std::shared_ptr< Nonlinear_Solver > Nonlinear_Solver_Factory::create_nonlinear_solver( const enum NonlinearSolverType   aNonLinSolverType )
{
    std::shared_ptr< Nonlinear_Solver > tNonLinSys;

    switch( aNonLinSolverType )
    {
    case ( NonlinearSolverType::NEWTON_SOLVER ):
        tNonLinSys = std::make_shared< Newton_Solver >();
        break;
    default:
        MORIS_ASSERT( false, "No solver type specified" );
        break;
    }

    return tNonLinSys;
}

