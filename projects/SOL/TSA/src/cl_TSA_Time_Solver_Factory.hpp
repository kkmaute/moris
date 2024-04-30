/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Factory.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_
#define SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_

#include <memory>

#include "cl_TSA_Time_Solver_Enums.hpp"

#include "cl_Parameter_List.hpp"

namespace moris
{
class Solver_Interface;
    namespace tsa
    {
        class Time_Solver_Algorithm;
        class Time_Solver_Factory
        {
        private:

        protected:

        public:
            Time_Solver_Factory();

            ~Time_Solver_Factory();

            std::shared_ptr< Time_Solver_Algorithm > create_time_solver( const enum TimeSolverType   aTimeSolverType = TimeSolverType::MONOLITHIC );

            std::shared_ptr< Time_Solver_Algorithm > create_time_solver( const enum TimeSolverType   aTimeSolverType,
                                                                         const Parameter_List                                                           aParameterlist);
        };
    }
}

#endif /* SRC_DISTLINALG_CL_TSA_TIME_SOLVER_FACTORY_HPP_ */

