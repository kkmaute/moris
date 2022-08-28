/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_Triangle_Vertex.cpp
 *
 */

#include "cl_SDF_Triangle_Vertex.hpp"
#include "op_times.hpp"

namespace moris
{
    namespace sdf
    {
//-------------------------------------------------------------------------------

        Triangle_Vertex::Triangle_Vertex(
                const moris_index        aIndex,
                const Matrix< DDRMat > & aNodeCoords ) :
                            mIndex( aIndex ),
                            mNodeCoords( aNodeCoords ),
                            mOriginalNodeCoords( aNodeCoords )
        {
        }

//-------------------------------------------------------------------------------

        void
        Triangle_Vertex::rotate_node_coords( const Matrix< F33RMat > & aRotationMatrix )
        {
            mNodeCoords = aRotationMatrix * mOriginalNodeCoords;
        }

//-------------------------------------------------------------------------------

        void
        Triangle_Vertex::reset_node_coords()
        {
            mNodeCoords = mOriginalNodeCoords;
        }

//-------------------------------------------------------------------------------

    } /* namespace sdf */
} /* namespace moris */

