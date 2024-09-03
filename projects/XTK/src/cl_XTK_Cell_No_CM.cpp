/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Cell_No_CM.cpp
 *
 */

#include "cl_XTK_Cell_No_CM.hpp"

#include <utility>
#include "cl_XTK_Background_Mesh.hpp"

namespace moris::xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------
    Cell_XTK_No_CM::Cell_XTK_No_CM( moris::moris_id aElementId,
            moris::moris_index                      aElementIndex,
            moris::moris_index                      aElementOwner,
            std::shared_ptr< mtk::Cell_Info >       aCellInfo,
            const Vector< mtk::Vertex* >&           aVertices )
            : Cell( aElementId, aElementIndex, aElementOwner, std::move( aCellInfo ) )
            , mCellVertices( aVertices )
    {
    }

    // ----------------------------------------------------------------------------------
    // Cell get functions
    // ----------------------------------------------------------------------------------

    Matrix< DDRMat >
    Cell_XTK_No_CM::get_vertex_coords() const
    {
        size_t           tNumVertices = this->get_number_of_vertices();
        Matrix< DDRMat > tVertexCoords;
        for ( size_t i = 0; i < tNumVertices; i++ )
        {
            Matrix< DDRMat > tVertCoord = mCellVertices( i )->get_coords();

            if ( i == 0 )
            {
                tVertexCoords.resize( tNumVertices, tVertCoord.numel() );
            }

            tVertexCoords.set_row( i, tVertCoord );
        }
        return tVertexCoords;
    }
}    // namespace moris::xtk
