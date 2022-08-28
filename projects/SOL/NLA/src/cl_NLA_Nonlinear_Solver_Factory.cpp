/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver_Factory.cpp
 *
 */

#include "cl_NLA_Nonlinear_Solver_Factory.hpp"
#include "cl_DLA_Solver_Interface.hpp"

#include "cl_NLA_Newton_Solver.hpp"
#include "cl_NLA_NLBGS.hpp"
#include "cl_NLA_Nonlinear_Algorithm.hpp"
//#include "cl_NLA_Arc_Length.hpp"

using namespace moris;
using namespace NLA;

Nonlinear_Solver_Factory::Nonlinear_Solver_Factory()
{}

Nonlinear_Solver_Factory::~Nonlinear_Solver_Factory()
{}

std::shared_ptr< Nonlinear_Algorithm > Nonlinear_Solver_Factory::create_nonlinear_solver( const enum NonlinearSolverType aNonLinSolverType )
{
    std::shared_ptr< Nonlinear_Algorithm > tNonLinSys = nullptr;

    switch( aNonLinSolverType )
    {
    case ( NonlinearSolverType::NEWTON_SOLVER ):
        tNonLinSys = std::make_shared< Newton_Solver >();
        break;
    case ( NonlinearSolverType::NLBGS_SOLVER ):
        tNonLinSys = std::make_shared< NonLinBlockGaussSeidel >();
        break;
//    case ( NonlinearSolverType::ARC_LENGTH_SOLVER ):
//        tNonLinSys = std::make_shared< Arc_Length_Solver >();
//        break;
    default:
        MORIS_ERROR( false, "No solver type specified" );
        break;
    }
    return tNonLinSys;
}

    std::shared_ptr< Nonlinear_Algorithm > Nonlinear_Solver_Factory::create_nonlinear_solver( const enum NonlinearSolverType aNonLinSolverType,
                                                                                              const ParameterList            aParameterlist )
    {
        std::shared_ptr< Nonlinear_Algorithm > tNonLinSys = nullptr;

        switch( aNonLinSolverType )
        {
        case ( NonlinearSolverType::NEWTON_SOLVER ):
            tNonLinSys = std::make_shared< Newton_Solver >( aParameterlist );
            break;
        case ( NonlinearSolverType::NLBGS_SOLVER ):
            tNonLinSys = std::make_shared< NonLinBlockGaussSeidel >( aParameterlist );
            break;
//        case ( NonlinearSolverType::ARC_LENGTH_SOLVER ):
//            tNonLinSys = std::make_shared< Arc_Length_Solver >();
//            break;
        default:
            MORIS_ERROR( false, "No solver type specified" );
            break;
        }

//    tNonLinSys->get_my_nonlin_solver()->set_nonlin_solver_type( aNonLinSolverType );

    return tNonLinSys;
}

