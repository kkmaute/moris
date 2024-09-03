/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Relaxation.hpp
 *
 */
#ifndef SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_
#define SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_

#include "cl_Parameter_List.hpp"

#include "cl_SOL_Enums.hpp"

namespace moris::NLA
{
    class Nonlinear_Solver;

    class Solver_Relaxation
    {
      private:
        /// relaxation strategy
        sol::SolverRelaxationType mRelaxationStrategy;

        /// basic relaxation parameter; use depends on strategy
        real mRelaxation = MORIS_REAL_MAX;

        /// damping factor for adaptation of relaxation parameter
        real mBetaDamping = MORIS_REAL_MAX;

        /// number of trials in adaptive search strategies
        uint mNumTrials = 0;

        /// parameters for adaptive strategies
        real mAlphak    = -1.0;
        real mBetak     = -1.0;
        real mResk      = -1.0;
        real mRelaxHist = MORIS_REAL_MAX;

      public:
        Solver_Relaxation( Parameter_List& aParameterListNonlinearSolver );

        ~Solver_Relaxation(){};

        /*
         *  evaluates the relaxation parameter
         */
        bool eval(
                sint              aIter,
                Nonlinear_Solver* aNonLinSolverManager,
                real&             aRelaxValue );
    };
}    // namespace moris::NLA

#endif /* SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_ */
