/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Amesos.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_
#define SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_

// TPL header files
#include "Epetra_ConfigDefs.h"
#include "Amesos_ConfigDefs.h"

#include "cl_DLA_Linear_Solver_Algorithm.hpp"

#include "Amesos.h"
#include "Amesos_BaseSolver.h"

namespace moris
{
    namespace dla
    {
        class Linear_Solver_Amesos : public Linear_Solver_Algorithm
        {
          private:
            Amesos_BaseSolver* mAmesosSolver = nullptr;

            Linear_Problem* mLinearSystem = nullptr;

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
    }    // namespace dla
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_LINEAR_SOLVER_AMESOS_HPP_ */
