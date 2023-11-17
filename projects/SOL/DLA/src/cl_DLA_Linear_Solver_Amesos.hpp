/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Amesos.hpp
 *
 */

#pragma once

// TPL header files
#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"

#include "cl_DLA_Linear_Solver_Algorithm_Trilinos.hpp"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

namespace moris::dla
{
    class Linear_Solver_Amesos : public Linear_Solver_Algorithm_Trilinos
    {
      private:
        Amesos_BaseSolver* mAmesosSolver = nullptr;

        Epetra_LinearProblem mEpetraProblem;

        bool mIsPastFirstSolve;

      protected:

      public:
        Linear_Solver_Amesos();

        Linear_Solver_Amesos( const moris::ParameterList aParameterlist );

        Linear_Solver_Amesos( Linear_Problem* aLinearSystem );

        ~Linear_Solver_Amesos();

        void set_solver_parameters();

        // int SetSystemMatrix ( bool aUseTranspose );

        moris::sint solve_linear_system();

        moris::sint solve_linear_system( Linear_Problem* aLinearSystem, const moris::sint aIter );

        void set_solver_internal_parameters();
    };
}    // namespace moris::dla
