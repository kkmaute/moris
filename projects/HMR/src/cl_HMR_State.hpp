/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_State.hpp
 *
 */

#pragma once

namespace moris::hmr
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
} /* namespace moris */
