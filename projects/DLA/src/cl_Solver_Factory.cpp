/*
 * cl_Solver_Factory.cpp
 *
 *  Created on: Apr 10, 2018
 *      Author: schmidt
 */
#include "cl_Solver_Factory.hpp"
#include "cl_Solver_Interface.hpp"

#include "cl_Linear_Solver_Trilinos.hpp"
#include "cl_Linear_Solver_Aztec.hpp"
#include "cl_Linear_Solver_Amesos.hpp"
#include "cl_Linear_Solver_Amesos2.hpp"
#include "cl_Linear_Solver_PETSc.hpp"

using namespace moris;

Solver_Factory::Solver_Factory()
{}

Solver_Factory::~Solver_Factory()
{}

std::shared_ptr< Linear_Solver > Solver_Factory::create_solver(       moris::Solver_Interface * aSolverInterface,
                                                                const enum SolverType           aSolverType )
{
    std::shared_ptr< Linear_Solver > tLinSys;

    switch( aSolverType )
    {
    case ( SolverType::TRILINOSTEST ):
        tLinSys = std::make_shared< Linear_Solver_Trilinos >( aSolverInterface );
        break;
    case ( SolverType::AZTEC_IMPL ):
        tLinSys = std::make_shared< Linear_Solver_Aztec >( aSolverInterface );
        break;
    case ( SolverType::AMESOS_IMPL ):
        tLinSys = std::make_shared< Linear_Solver_Amesos >( aSolverInterface );
        break;
    case ( SolverType::AMESOS2_IMPL ):
        tLinSys = std::make_shared< Linear_Solver_Amesos2 >( aSolverInterface );
        break;
    case ( SolverType::PETSC):
        tLinSys = std::make_shared< Linear_Solver_PETSc >( aSolverInterface );
        break;
    default:
        MORIS_ASSERT( false, "No solver type specified" );
        break;
    }


    return tLinSys;
}


