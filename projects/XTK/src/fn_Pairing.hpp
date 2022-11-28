/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_Pairing.hpp
 *
 */

#ifndef SRC_TOOLS_FN_PAIRING_HPP_
#define SRC_TOOLS_FN_PAIRING_HPP_

#include "cl_Matrix.hpp"

namespace xtk
{
    template< typename Integer >
    inline static Integer
    cantor_pairing( Integer const aPt1,
            Integer const         aPt2 )
    {
        Integer tID = ( ( aPt1 * aPt1 + 3 * aPt1 + 2 * aPt1 * aPt2 + aPt2 + aPt2 * aPt2 ) / 2 );
        return tID;
    }
}    // namespace xtk

#endif /* SRC_TOOLS_FN_PAIRING_HPP_ */
