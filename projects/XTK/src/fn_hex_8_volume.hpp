/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_hex_8_volume.hpp
 *
 */

#ifndef PROJECTS_XTK_SRC_TOOLS_FN_HEX_8_VOLUME_HPP_
#define PROJECTS_XTK_SRC_TOOLS_FN_HEX_8_VOLUME_HPP_

#include "cl_Matrix.hpp"

namespace xtk
{
    /*
     * Compute the volume of a hexahedron
     */

    template< typename Real_Matrix >
    inline typename moris::Matrix< Real_Matrix >::Data_Type
    compute_hex_8_volume(
            moris::Matrix< Real_Matrix > const &     aNodeCoordinates,
            moris::Matrix< moris::IndexMat > const & aNodeToElement )
    {
        typedef typename moris::Matrix< Real_Matrix >::Data_Type Real;

        // Choose nodes 0 and 7 to get lengths
        Real tLx = std::abs( aNodeCoordinates( aNodeToElement( 0 ), 0 ) - aNodeCoordinates( aNodeToElement( 6 ), 0 ) );
        Real tLy = std::abs( aNodeCoordinates( aNodeToElement( 0 ), 1 ) - aNodeCoordinates( aNodeToElement( 6 ), 1 ) );
        Real tLz = std::abs( aNodeCoordinates( aNodeToElement( 0 ), 2 ) - aNodeCoordinates( aNodeToElement( 6 ), 2 ) );

        return tLx * tLy * tLz;
    }
}    // namespace xtk

#endif /* PROJECTS_XTK_SRC_TOOLS_FN_HEX_8_VOLUME_HPP_ */
