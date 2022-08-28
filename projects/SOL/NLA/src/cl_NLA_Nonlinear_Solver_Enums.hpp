/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_NLA_Nonlinear_Solver_Enums.hpp
 *
 */

#ifndef SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_ENUMS_HPP_
#define SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_ENUMS_HPP_

namespace moris
{
    namespace NLA
    {
        enum class NonlinearSolverType
        {
            NEWTON_SOLVER,   //< Wrapper around Aztec Solver
            NLBGS_SOLVER,   //< Wrapper around Aztec Solver
            ARC_LENGTH_SOLVER,
            END_ENUM
        };
    }
}

#endif /* SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_ENUMS_HPP_ */

