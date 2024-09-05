/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_State.hpp
 *
 */

#ifndef PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_
#define PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    enum class State
    {
        PRINT_USAGE,
        PRINT_HELP,
        PRINT_VERSION,
        CALCULATE_RAYCAST,
        CALCULATE_RAYCAST_AND_SDF
    };

    //-------------------------------------------------------------------------------
    }

#endif /* PROJECTS_GEN_SDF_SRC_CL_SDF_STATE_HPP_ */

