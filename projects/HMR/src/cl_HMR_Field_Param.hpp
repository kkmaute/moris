/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Field_Param.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"

namespace moris::hmr
{
// -----------------------------------------------------------------------------

    struct Field_Param
    {
        std::string mLabel;
        moris_id    mID = gNoID;
        uint        mInputLagrangeOrder = 0;
        uint        mInputBSplineOrder = 0;
        uint        mOutputBSplineOrder = 0;
        std::string mSource;
        std::string mTarget;
        real        mL2alpha = 0.0;

        Field_Param(){};

        ~Field_Param(){};
    };

// -----------------------------------------------------------------------------
}
