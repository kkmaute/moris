/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_ML.hpp
 *
 */

#pragma once

#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris::dla
{
    class Linear_Problem;
    class Linear_Solver_ML : public Linear_Solver_Algorithm_Trilinos
    {

      private:
        // -----------------------------------------------------------------------------------

        bool build_external_preconditioner( const moris::sint &aIter = 1 );

        // -----------------------------------------------------------------------------------

      public:

        //-----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Linear_Solver_ML object
         * 
         * @param aParameterList
         */
        Linear_Solver_ML( const Parameter_List& aParameterList = prm::create_linear_algorithm_parameter_list( sol::SolverType::ML ) );

        // -----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Linear_Solver_ML object
         * 
         */
        ~Linear_Solver_ML() override = default;

        // -----------------------------------------------------------------------------------

        /**
         * @brief solve linear system by iterating
         * 
         * @param aLinearSystem 
         * @param aIter 
         * @return moris::sint 
         */
        moris::sint solve_linear_system(
                Linear_Problem   *aLinearSystem,
                const moris::sint aIter ) override;

        //-----------------------------------------------------------------------------------

        /**
         * @brief solve linear system not implemented yet
         * 
         * @return moris::sint 
         */
        moris::sint solve_linear_system() override;
    };
}    // namespace moris::dla
