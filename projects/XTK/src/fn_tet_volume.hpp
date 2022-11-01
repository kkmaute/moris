/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_tet_volume.hpp
 *
 */

#ifndef SRC_TOOLS_FN_TET_VOLUME_HPP_
#define SRC_TOOLS_FN_TET_VOLUME_HPP_

#include "cl_Matrix.hpp"
#include "fn_det.hpp"

namespace xtk
{
    inline moris::real
    vol_tetrahedron( moris::Matrix< moris::DDRMat >& aCoord )
    {
        // explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
        moris::Matrix< moris::DDRMat > J( 4, 4, 1.0 );    // 4 Coordinates are needed for the Jacobian-matrix.
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

        typename moris::Matrix< moris::DDRMat >::Data_Type volume = det( J ) / 6.0;    // Volume = 1/6*det(J)

        return volume;
    }

    //------------------------------------------------------------------------------

    inline moris::real
    vol_tetrahedron( moris::Matrix< moris::DDRMat > const & aCoord,
            moris::Matrix< moris::IndexMat > const &        aNodeToElement )
    {
        // explanation: www.colorado.edu/engineering/Aerospace/CAS/courses.d/AFEM.d/AFEM.Ch09.d/AFEM.Ch09.pdf
        moris::Matrix< moris::DDRMat > J( 4, 4, 1.0 );    // 4 Coordinates are needed for the Jacobian-matrix.
        // The order of the coordinates is irrelevant.

        J( 1, 0 ) = aCoord( aNodeToElement( 0 ), 0 );
        J( 1, 1 ) = aCoord( aNodeToElement( 1 ), 0 );
        J( 1, 2 ) = aCoord( aNodeToElement( 2 ), 0 );
        J( 1, 3 ) = aCoord( aNodeToElement( 3 ), 0 );
        J( 2, 0 ) = aCoord( aNodeToElement( 0 ), 1 );
        J( 2, 1 ) = aCoord( aNodeToElement( 1 ), 1 );
        J( 2, 2 ) = aCoord( aNodeToElement( 2 ), 1 );
        J( 2, 3 ) = aCoord( aNodeToElement( 3 ), 1 );
        J( 3, 0 ) = aCoord( aNodeToElement( 0 ), 2 );
        J( 3, 1 ) = aCoord( aNodeToElement( 1 ), 2 );
        J( 3, 2 ) = aCoord( aNodeToElement( 2 ), 2 );
        J( 3, 3 ) = aCoord( aNodeToElement( 3 ), 2 );

        typename moris::Matrix< moris::DDRMat >::Data_Type volume = det( J ) / 6.0;    // Volume = 1/6*det(J)

        return volume;
    }

}    // namespace xtk

#endif /* SRC_TOOLS_FN_TET_VOLUME_HPP_ */
