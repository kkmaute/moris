/*
 * cl_NLA_Solver_Relaxation.hpp
 *
 *  Created on: Sep 21, 2018
 *      Author: schmidt
 */
#ifndef SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_
#define SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_

#include "cl_Param_List.hpp"

#include "cl_SOL_Enums.hpp"

namespace moris
{
    namespace NLA
    {
        class Solver_Relaxation
        {
            private:

                /// relaxation strategy
                sol::SolverRelaxationType mRelaxationStratey;

                /// basic relaxation parameter; use depends on strategy
                real mRelaxation = MORIS_REAL_MAX;

                /// damping factor for adaptation of relaxation parameter
                real mBetaDamping =  MORIS_REAL_MAX;

                /// number of trials in adaptive search strategies
                uint mNumTrials = 0;

                /// parameters for adaptive strategies
                real mAlphak      = -1.0;
                real mBetak       = -1.0;
                real mResk        = -1.0;
                real mRelaxHist   = MORIS_REAL_MAX;

            public:
                Solver_Relaxation( ParameterList & aParameterListNonlinearSolver);

                ~Solver_Relaxation(){};

                /*
                 *  evaluates the relaxation parameter
                 */
                bool eval(
                        const real & aRefNorm,
                        const real & aResNorm,
                        real       & aRelaxValue );
        };
    }
}

#endif /* SRC_FEM_CL_NLA_SOLVER_RELAXATION_HPP_ */
