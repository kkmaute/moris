/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Solver_Factory.hpp
 *
 */

#pragma once

#include <memory>

#include "cl_SOL_Enums.hpp"

#include "cl_Parameter_List.hpp"

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
        class Preconditioner;

        class Solver_Factory
        {
          private:

          protected:

          public:
            Solver_Factory();

            ~Solver_Factory();

            //------------------------------------------------------------------------------

            /**
             * @brief
             *
             * @param aPreconditionerType
             * @param aParameterlist
             * @return Preconditioner_Trilinos*
             */
            Preconditioner*
            create_preconditioner( const Parameter_List& aParameterlist );

            /**
             * Creates a solver of the given type and with given parameters.
             *
             * @param aSolverType Solver type
             * @param aParameterlist Solve parameter list
             * @return Solver shared pointer
             */
            std::shared_ptr< Linear_Solver_Algorithm > create_solver(
                    const Parameter_List& aParameterlist );

            //------------------------------------------------------------------------------

            Linear_Problem* create_linear_system(
                    moris::Solver_Interface* aSolverInterface,
                    sol::MapType             aLinSysType               = sol::MapType::Epetra,
                    bool                     aNotCreatedByNonLinSolver = false );
            
            //------------------------------------------------------------------------------
            
            Linear_Problem* create_linear_system(
                    moris::Solver_Interface* aSolverInterface,
                    sol::SOL_Warehouse*      aSolverWarehouse,
                    sol::Dist_Map*           aMap,
                    sol::Dist_Map*           aFullMap,
                    sol::MapType             aLinSysType               = sol::MapType::Epetra,
                    bool                     aNotCreatedByNonLinSolver = false );
        };
    }    // namespace dla
}    // namespace moris

