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

#include "assert.hpp"
#include "cl_Map.hpp"
#include "fn_enum_macros.hpp"
#include "cl_Vector.hpp"


namespace moris
{
    namespace NLA
    {
        ENUM_MACRO( NonlinearSolverType,
                NEWTON_SOLVER,
                NLBGS_SOLVER,
                ARC_LENGTH_SOLVER,
                END_ENUM )
    }
}    // namespace moris

#endif /* SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_ENUMS_HPP_ */
