/*
 * cl_NLA_Nonlinear_Solver_Factory.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#include <memory>

#include "cl_NLA_Nonlinear_Solver.hpp"

#ifndef SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_

namespace moris
{
class Solver_Interface;
    namespace NLA
    {
        class Nonlinear_Solver_Factory
        {
        private:

        protected:

        public:
            Nonlinear_Solver_Factory();

            ~Nonlinear_Solver_Factory();

            std::shared_ptr< Nonlinear_Solver > create_nonlinear_solver( Solver_Interface               * aSolverInput,
                                                                         const enum NonlinearSolverType   aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );

            std::shared_ptr< Nonlinear_Solver > create_nonlinear_solver( const enum NonlinearSolverType   aNonLinSolverType = NonlinearSolverType::NEWTON_SOLVER );
        };
    }
}

#endif /* SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_FACTORY_HPP_ */
