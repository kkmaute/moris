/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Facet_Vertex.cpp
 *
 */

#include "cl_SDF_Facet_Vertex.hpp"
#include "op_times.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    Facet_Vertex::Facet_Vertex(
            const moris_index       aIndex,
            const Matrix< DDRMat > &aNodeCoords )
            : mIndex( aIndex )
            , mNodeCoords( aNodeCoords )
            , mOriginalNodeCoords( aNodeCoords )
    {
    }

    //-------------------------------------------------------------------------------

    void
    Facet_Vertex::rotate_node_coords( const Matrix< DDRMat > &aRotationMatrix )
    {
        mNodeCoords    = aRotationMatrix * mNodeCoords;
        mIsTransformed = true;
    }

    //-------------------------------------------------------------------------------

    void
    Facet_Vertex::reset_node_coords()
    {
        mNodeCoords = mOriginalNodeCoords;
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf
