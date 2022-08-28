/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver_Factory.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_

#include <memory>

#include "cl_NLA_Nonlinear_Solver_Enums.hpp"

#include "cl_Param_List.hpp"

namespace moris
{
class Solver_Interface;
    namespace NLA
    {
        class Nonlinear_Algorithm;
        class Nonlinear_Solver_Factory
        {
        private:

        protected:

        public:
            Nonlinear_Solver_Factory();

            ~Nonlinear_Solver_Factory();

            std::shared_ptr< Nonlinear_Algorithm > create_nonlinear_solver( const enum NonlinearSolverType   aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

            std::shared_ptr< Nonlinear_Algorithm > create_nonlinear_solver( const enum NonlinearSolverType   aNonLinSolverType,
                                                                            const ParameterList              aParameterlist );
        };
    }
}

#endif /* SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_ */

