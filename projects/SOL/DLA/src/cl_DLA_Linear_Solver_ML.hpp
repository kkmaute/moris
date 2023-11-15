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

#include "cl_Param_List.hpp"    //CNT/src

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
        // -----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Linear_Solver_ML object
         * 
         */
        Linear_Solver_ML(){};

        //-----------------------------------------------------------------------------------

        /**
         * @brief Construct a new Linear_Solver_ML object
         * 
         * @param aParameterlist 
         */
        Linear_Solver_ML( const moris::ParameterList aParameterlist );

        // -----------------------------------------------------------------------------------

        /**
         * @brief Destroy the Linear_Solver_ML object
         * 
         */
        virtual ~Linear_Solver_ML();


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
                const moris::sint aIter );

        //-----------------------------------------------------------------------------------

        /**
         * @brief solve linear system not implemented yet
         * 
         * @return moris::sint 
         */
        moris::sint solve_linear_system();
    };
}    // namespace moris::dla
