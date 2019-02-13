/*
 * cl_TSA_Time_Solver_Enums.hpp
 *
 *  Created on: Jul 2, 2018
 *      Author: schmidt
 */
#ifndef SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_
#define SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_

namespace moris
{
    namespace tsa
    {
        enum class TimeSolverType
        {
            MONOLITHIC,   //< Wrapper around Aztec Solver
            STAGGERED,   //< Wrapper around Aztec Solver
            END_ENUM
        };
    }
}

#endif /* SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_ */
