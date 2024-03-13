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

namespace moris::NLA
{
    std::shared_ptr< Nonlinear_Algorithm > Nonlinear_Solver_Factory::create_nonlinear_solver( const Parameter_List& aParameterList )
    {
        std::shared_ptr< Nonlinear_Algorithm > tNonLinSys = nullptr;

        switch ( static_cast< NonlinearSolverType >( aParameterList.get< uint >( "NLA_Solver_Implementation" ) ) )
        {
            case ( NonlinearSolverType::NEWTON_SOLVER ):
                tNonLinSys = std::make_shared< Newton_Solver >( aParameterList );
                break;
            case ( NonlinearSolverType::NLBGS_SOLVER ):
                tNonLinSys = std::make_shared< NonLinBlockGaussSeidel >( aParameterList );
                break;
        //        case ( NonlinearSolverType::ARC_LENGTH_SOLVER ):
        //            tNonLinSys = std::make_shared< Arc_Length_Solver >();
        //            break;
            default:
                MORIS_ERROR( false, "No solver type specified" );
                break;
        }

        return tNonLinSys;
    }
}
