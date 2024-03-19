/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Load_Control.hpp
 *
 */
#ifndef SRC_FEM_CL_NLA_SOLVER_LOAD_CONTROL_HPP_
#define SRC_FEM_CL_NLA_SOLVER_LOAD_CONTROL_HPP_

#include "cl_Param_List.hpp"

#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace NLA
    {
        class Nonlinear_Solver;

        class Solver_Load_Control
        {
          private:
            /// relaxation strategy
            sol::SolverLoadControlType mLoadControlStrategy;

            /// load step counter
            uint mLoadStepCounter;

            /// maximum number of load steps
            uint mNumLoadSteps;

            /// initial (constant) load factor
            real mInitialLoadFactor;

            /// required relative residual drop to increase load factor
            real mRelativeResidualDropThreshold;

            /// parameters for exponential strategy
            real mExponent;

            // check if the load stepping requirement (e.g. residual drop) is met
            bool check_load_step_requirement( Nonlinear_Solver* aNonLinSolverManager );

          public:
            Solver_Load_Control( ParameterList& aParameterListNonlinearSolver );

            ~Solver_Load_Control(){};

            /*
             *  evaluates the relaxation parameter
             */
            real get_initial_load_factor();

            /*
             *  evaluates the relaxation parameter
             */
            void eval(
                    sint              aIter,
                    Nonlinear_Solver* aNonLinSolverManager,
                    real&             aLoadFactor );
        };
    }    // namespace NLA
}    // namespace moris

#endif /* SRC_FEM_CL_NLA_SOLVER_LOAD_CONTROL_HPP_ */
