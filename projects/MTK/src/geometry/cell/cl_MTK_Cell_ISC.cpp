/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_ISC.cpp
 *
 */

#include "cl_MTK_Cell_ISC.hpp"

namespace moris
{

    namespace mtk
    {
        // ----------------------------------------------------------------------------------
        // Constructor/Deconstructor Source code
        // ----------------------------------------------------------------------------------
        Cell_ISC::Cell_ISC(moris::moris_id                       aCellId,
                moris::moris_index              aCellIndex,
                moris::moris_id                 aCellOwner,
                std::shared_ptr<mtk::Cell_Info> aCellInfo,
                moris::Vector< mtk::Vertex* >     aVertices)
                                           : Cell(aCellId,aCellIndex,aCellOwner, aCellInfo),
                                             mCellVertices(aVertices)
                                             {}

        // ----------------------------------------------------------------------------------
        // Cell get functions
        // ----------------------------------------------------------------------------------
        Matrix< DDRMat >
        Cell_ISC::get_vertex_coords() const
        {
            size_t tNumVertices = this->get_number_of_vertices();
            Matrix< DDRMat > tVertexCoords;
            for(size_t i = 0; i<tNumVertices; i++)
            {
                Matrix<DDRMat> tVertCoord = mCellVertices(i)->get_coords();

                if(i == 0 )
                {
                    tVertexCoords.resize(tNumVertices,tVertCoord.numel());
                }

                tVertexCoords.set_row(i,tVertCoord);
            }
            return tVertexCoords;
        }
    }
}

