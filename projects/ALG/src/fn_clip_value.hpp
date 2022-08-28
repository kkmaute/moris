/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_clip_value.hpp
 *
 */

#ifndef SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_
#define SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_

#include <iostream>

// MORIS header files.
#include "typedefs.hpp"    // COR/src
#include "linalg_typedefs.hpp"

namespace moris
{
    /**
     * clip_value function
     *
     * clips value based on threshold, i.e., small value close to zero maintaining sign
     *
     *    if abs(value) < threshold && value >= 0 :  threshold
     *    if abs(value) < threshold && value <  0 : -threshold
     *    otherwise:                                 aValue
     *
     * @param[in] aValue     value to be clipped
     * @param[in] aThreshold threshold, default: MORIS_REAL_EPS
     *
     * @return clipped value
     */
    inline real
    clip_value(
            const real& aValue,
            const real& aThreshold = MORIS_REAL_EPS )
    {
        if ( std::abs( aValue ) < aThreshold )
        {
            return aValue < 0.0 ? -aThreshold : aThreshold;
        }

        return aValue;
    }
}    // namespace moris

#endif /* SRC_ALGORITHMS_FN_CLIP_VALUE_HPP_ */

