/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Solver_Nonconformal_Remapping.cpp
 *
 */

#pragma once

#include "cl_SOL_Enums.hpp"
#include "moris_typedefs.hpp"

namespace moris::NLA
{
    class Nonlinear_Solver;
    class Solver_Nonconformal_Remapping
    {
      public:
        Solver_Nonconformal_Remapping( Parameter_List& aParameterListNonlinearSolver );

        bool requires_remapping( uint aIter, Nonlinear_Solver* aNonLinSolverManager, real& aLoadFactor );

      private:
        sol::SolverRaytracingStrategy mStrategy;
        uint                          mLoadStepFrequency{};
        uint                          mIterationFrequency{};
        real                          mResidualChangeTolerance{};
        real                          mPreviousResidual{ -1.0 };
        real                          mPreviousLoadFactor;
        uint                          mLoadStepCounter;

        bool check_every_nth_load_step( real& aLoadFactor );

        static bool check_every_nth_iteration( uint aIter, uint aFrequency );

        bool check_relative_residual_change( Nonlinear_Solver* aNonLinSolverManager );
    };

}    // namespace moris::NLA
