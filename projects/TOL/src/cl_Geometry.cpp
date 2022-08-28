/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Geometry.cpp
 *
 */

#include "cl_Geometry.hpp"
#include "algorithms.hpp"
#include "fn_det.hpp"

namespace moris
{
    real
    Geometry::vol_tetrahedron( Matrix< DDRMat >& aCoord )
    {
        // explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
        Matrix< DDRMat > J( 4, 4, 1.0 );    // 4 Coordinates are needed for the Jacobinan-matrix.
        // The order of the coordinates is irrelevant.

        J( 1, 0 ) = aCoord( 0, 0 );
        J( 1, 1 ) = aCoord( 1, 0 );
        J( 1, 2 ) = aCoord( 2, 0 );
        J( 1, 3 ) = aCoord( 3, 0 );
        J( 2, 0 ) = aCoord( 0, 1 );
        J( 2, 1 ) = aCoord( 1, 1 );
        J( 2, 2 ) = aCoord( 2, 1 );
        J( 2, 3 ) = aCoord( 3, 1 );
        J( 3, 0 ) = aCoord( 0, 2 );
        J( 3, 1 ) = aCoord( 1, 2 );
        J( 3, 2 ) = aCoord( 2, 2 );
        J( 3, 3 ) = aCoord( 3, 2 );

        real volume = std::abs( moris::det( J ) / 6.0 );    // Volume = 1/6*det(J)

        return volume;
    }
}    // namespace moris

