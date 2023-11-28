/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Aztec.hpp
 *
 */

#pragma once

#include "AztecOO.h"

#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"

#include "cl_Param_List.hpp"    //CNT/src

namespace moris
{
    namespace dla
    {
        class Linear_Problem;
        class Linear_Solver_Aztec : public Linear_Solver_Algorithm_Trilinos
        {
          private:

            AztecOO *mAztecSolver = nullptr;

            Epetra_LinearProblem mEpetraProblem;

          private:
            // -----------------------------------------------------------------------------------

            void set_solver_internal_parameters();

            // -----------------------------------------------------------------------------------

            void set_solver_parameters();

            // -----------------------------------------------------------------------------------

            bool build_external_preconditioner( const moris::sint &aIter = 1 );

            // -----------------------------------------------------------------------------------

          public:
            // -----------------------------------------------------------------------------------

            Linear_Solver_Aztec();

            // -----------------------------------------------------------------------------------

            Linear_Solver_Aztec( const moris::ParameterList aParameterlist );

            // -----------------------------------------------------------------------------------

            Linear_Solver_Aztec( Linear_Problem *aLinearSystem );

            // -----------------------------------------------------------------------------------

            ~Linear_Solver_Aztec();

            // -----------------------------------------------------------------------------------

            moris::sint solve_linear_system();

            // -----------------------------------------------------------------------------------

            moris::sint solve_linear_system(
                    Linear_Problem   *aLinearSystem,
                    const moris::sint aIter );

            // -----------------------------------------------------------------------------------
        };
    }    // namespace dla
}    // namespace moris

