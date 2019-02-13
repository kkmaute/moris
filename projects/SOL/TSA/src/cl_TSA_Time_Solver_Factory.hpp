/*
 * cl_TSA_Time_Solver_Factory.hpp
 *
 *  Created on: Mar 28, 2018
 *      Author: schmidt
 */
#include <memory>

#include "cl_TSA_Time_Solver.hpp"
#include "cl_TSA_Time_Solver_Enums.hpp"

#ifndef SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_

namespace moris
{
class Solver_Interface;
    namespace tsa
    {
        class Time_Solver_Factory
        {
        private:

        protected:

        public:
            Time_Solver_Factory();

            ~Time_Solver_Factory();

            Time_Solver * create_time_solver( const enum TimeSolverType   aTimeSolverType = TimeSolverType::MONOLITHIC );
        };
    }
}

#endif /* SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_ */
