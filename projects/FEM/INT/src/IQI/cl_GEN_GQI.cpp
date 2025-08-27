/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_GQI.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_GEN_GQI.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    GQI::GQI() {}

    //------------------------------------------------------------------------------

    void
    GQI::initialize()
    {
        if ( !mIsInitialized )
        {
            // set initialize flag to true
            mIsInitialized = true;
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
