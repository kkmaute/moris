/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Solver_Factory.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_

#include <memory>

#include "cl_SOL_Enums.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
    class Solver_Interface;

    namespace sol
    {
        class SOL_Warehouse;
        class Dist_Map;
    }    // namespace sol
    namespace dla
    {
        class Linear_Solver_Algorithm;
        class Linear_Problem;

        class Solver_Factory
        {
          private:

          protected:

          public:
            Solver_Factory();

            ~Solver_Factory();

            std::shared_ptr< Linear_Solver_Algorithm > create_solver( const enum sol::SolverType aSolverType = sol::SolverType::AZTEC_IMPL );

            std::shared_ptr< Linear_Solver_Algorithm > create_solver(
                    const enum sol::SolverType aSolverType,
                    const ParameterList        aParameterlist );

            Linear_Problem* create_linear_system(
                    moris::Solver_Interface* aSolverInterface,
                    const enum sol::MapType  aLinSysType               = sol::MapType::Epetra,
                    const bool               aNotCreatedByNonLinSolver = false );

            Linear_Problem* create_linear_system(
                    moris::Solver_Interface* aSolverInterface,
                    sol::SOL_Warehouse*      aSolverWarehouse,
                    sol::Dist_Map*           aMap,
                    sol::Dist_Map*           aFullMap,
                    const enum sol::MapType  aLinSysType               = sol::MapType::Epetra,
                    const bool               aNotCreatedByNonLinSolver = false );
        };
    }    // namespace dla
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_DLA_SOLVER_FACTORY_HPP_ */
