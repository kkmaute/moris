/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_State.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_

namespace moris
{
    namespace hmr
    {
// -----------------------------------------------------------------------------

        enum class State
        {
            PRINT_USAGE,
            PRINT_HELP,
            PRINT_VERSION,
            INITIALIZE_MESH,
            REFINE_MESH,
            MAP_FIELDS
        };
// -----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_STATE_HPP_ */

