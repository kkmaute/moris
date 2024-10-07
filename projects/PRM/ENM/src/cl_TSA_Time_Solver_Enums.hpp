/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_TSA_Time_Solver_Enums.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_
#define SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_

namespace moris::tsa
{
    enum class TimeSolverType
    {
        MONOLITHIC,
        STAGGERED,
        END_ENUM
    };
}    // namespace moris::tsa

#endif /* SRC_DISTLINALG_CL_TSA_TIME_SOLVER_ENUMS_HPP_ */
