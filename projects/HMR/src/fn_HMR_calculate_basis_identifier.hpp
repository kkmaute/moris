/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_calculate_basis_identifier.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"

namespace moris::hmr
{
    /**
     * Calculates a unique basis identifier for a given IJK position and offsets.
     * Can be an index or an ID depending on the offset.
     *
     * @tparam N Number of dimensions of the IJK and offset arrays
     * @param aIJK IJK position
     * @param aDimensionOffset Offset array for each dimension
     * @return Unique basis identifier
     */
    template< uint N >
    luint calculate_basis_identifier(
            const luint* aIJK,
            const luint* aDimensionOffset )
    {
        luint tIdentifier = 0;
        for ( uint iDimension = 0; iDimension < N; iDimension++)
        {
            luint tOffsetTerm = aIJK[ iDimension ];
            for ( uint iPreviousDimension = 0; iPreviousDimension < iDimension; iPreviousDimension++ )
            {
                tOffsetTerm *= aDimensionOffset[ iPreviousDimension ];
            }
            tIdentifier += tOffsetTerm;
        }
        return tIdentifier;
    }
}
