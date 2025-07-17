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

#include "fn_enum_macros.hpp"


namespace moris::NLA
{
    ENUM_MACRO( NonlinearSolverType,
            NEWTON_SOLVER,
            NLBGS_SOLVER,
            ARC_LENGTH_SOLVER,
            END_ENUM )
}

#endif /* SRC_DISTLINALG_CL_NLA_NONLINEAR_SOLVER_ENUMS_HPP_ */
