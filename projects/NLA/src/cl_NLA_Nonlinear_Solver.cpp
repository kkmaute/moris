/*
 * cl_NLA_Nonlinear_Solver.cpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#include "cl_NLA_Nonlinear_Solver.hpp"

#include <ctime>

#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"

//#include "cl_DLA_Linear_Solver_Manager.hpp"
#include "cl_Vector.hpp"

#include "cl_Communication_Tools.hpp"

using namespace moris;
using namespace NLA;
using namespace dla;

//Nonlinear_Solver::Nonlinear_Solver()
//{
//}

Nonlinear_Solver::Nonlinear_Solver( Solver_Interface * aSolverInterface ) : mLinSolverManager()
{
    mNonlinearProblem = new Nonlinear_Problem( aSolverInterface );
}

void Nonlinear_Solver::set_linear_solver( std::shared_ptr< dla::Linear_Solver > aLinearSolver )
{
    mLinSolverManager.set_linear_solver( aLinearSolver );
}

void Nonlinear_Solver::set_linear_solver( const moris::uint aListEntry,
                                                std::shared_ptr< dla::Linear_Solver > aLinearSolver )
{
    mLinSolverManager.set_linear_solver( aListEntry, aLinearSolver );
}

void Nonlinear_Solver::set_nonlinear_problem( Nonlinear_Problem * aNonlinearProblem )
{
    mNonlinearProblem = aNonlinearProblem;
}


