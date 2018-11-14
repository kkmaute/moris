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
#include "cl_DLA_Linear_Solver_PETSc.hpp"

#include "cl_DLA_Linear_System_Trilinos.hpp"

using namespace moris;
using namespace dla;

Solver_Factory::Solver_Factory()
{}

Solver_Factory::~Solver_Factory()
{}

std::shared_ptr< Linear_Solver > Solver_Factory::create_solver( const enum SolverType aSolverType )
{
    std::shared_ptr< Linear_Solver > tLinSol;

    switch( aSolverType )
    {
    case ( SolverType::AZTEC_IMPL ):
        tLinSol = std::make_shared< Linear_Solver_Aztec >();
        break;
//    case ( SolverType::AMESOS_IMPL ):
//        tLinSol = std::make_shared< Linear_Solver_Amesos >( aSolverInterface );
//        break;
//    case ( SolverType::AMESOS2_IMPL ):
//        tLinSol = std::make_shared< Linear_Solver_Amesos2 >( aSolverInterface );
//        break;
//    case ( SolverType::PETSC):
//        tLinSol = std::make_shared< Linear_Solver_PETSc >( aSolverInterface );
//        break;
    default:
        MORIS_ASSERT( false, "No solver type specified" );
        break;
    }
    return tLinSol;
}

    std::shared_ptr< Linear_Problem > Solver_Factory::create_linear_system(       moris::Solver_Interface * aSolverInterface,
                                                                           const enum MapType              aLinSysType )
    {
        std::shared_ptr< Linear_Problem > tLinSys;

        switch( aLinSysType )
        {
        case ( MapType::Epetra ):
            tLinSys = std::make_shared< Linear_System_Trilinos >( aSolverInterface );
            break;
//        case ( MapType::Petsc):
//            tLinSys = std::make_shared< Linear_Solver_PETSc >( aSolverInterface );
//            break;
        default:
            MORIS_ASSERT( false, "No solver type specified" );
            break;
        }


    return tLinSys;
}


