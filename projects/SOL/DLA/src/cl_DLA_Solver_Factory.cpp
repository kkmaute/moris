/*
 * cl_DLA_Solver_Factory.cpp
 *
 *  Created on: Apr 10, 2018
 *      Author: schmidt
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

using namespace moris;
using namespace dla;

Solver_Factory::Solver_Factory()
{}

Solver_Factory::~Solver_Factory()
{}

std::shared_ptr< Linear_Solver_Algorithm > Solver_Factory::create_solver( const enum SolverType aSolverType )
{
    std::shared_ptr< Linear_Solver_Algorithm > tLinSol;

    switch( aSolverType )
    {
    case ( SolverType::AZTEC_IMPL ):
        tLinSol = std::make_shared< Linear_Solver_Aztec >();
        break;
    case ( SolverType::AMESOS_IMPL ):
        tLinSol = std::make_shared< Linear_Solver_Amesos >();
        break;
    case ( SolverType::BELOS_IMPL ):
        tLinSol = std::make_shared< Linear_Solver_Belos >();
        break;
    case ( SolverType::AMESOS2_IMPL ):
        tLinSol = std::make_shared< Linear_Solver_Amesos2 >();
        break;
    case ( SolverType::PETSC):
        tLinSol = std::make_shared< Linear_Solver_PETSc >(  );
        break;
    default:
        MORIS_ERROR( false, "No solver type specified" );
        break;
    }
    return tLinSol;
}

//-------------------------------------------------------------------------------------------------------------

Linear_Problem * Solver_Factory::create_linear_system(       moris::Solver_Interface * aSolverInterface,
		Dist_Map               * aMap,
		Dist_Map               * aFullMap,
                                                       const enum MapType              aLinSysType,
                                                       const bool                      aNotCreatedByNonLinSolver )
{
    Linear_Problem * tLinSys;

    switch( aLinSysType )
    {
    case ( MapType::Epetra ):
        tLinSys = new Linear_System_Trilinos( aSolverInterface, aMap, aFullMap );
        break;
    case ( MapType::Petsc):
        tLinSys = new Linear_System_PETSc( aSolverInterface, aMap, aFullMap, aNotCreatedByNonLinSolver );
        break;
    default:
        MORIS_ERROR( false, "Solver_Factory::create_linear_system: No solver type specified" );
        break;
    }

return tLinSys;
}

//--------------------------------------------------------------------------------------------------------------

Linear_Problem * Solver_Factory::create_linear_system(       moris::Solver_Interface * aSolverInterface,
                                                       const enum MapType              aLinSysType,
                                                       const bool                      aNotCreatedByNonLinSolver )
{
    Linear_Problem * tLinSys = nullptr;

    switch( aLinSysType )
    {
    case ( MapType::Epetra ):
        tLinSys = new Linear_System_Trilinos( aSolverInterface );
        break;
    case ( MapType::Petsc):
        tLinSys = new Linear_System_PETSc( aSolverInterface, aNotCreatedByNonLinSolver );
        break;
    default:
        MORIS_ERROR( false, "Solver_Factory::create_linear_system: No solver type specified" );
        break;
    }

return tLinSys;
}

