/*
 * cl_MTK_Cell_XTK_CM_Impl.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: doble
 */

#include "cl_XTK_Cell_No_CM.hpp"
#include "cl_XTK_Background_Mesh.hpp"

namespace xtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------
    Cell_XTK_No_CM::Cell_XTK_No_CM(moris::moris_id              aElementId,
                                   moris::moris_index          aElementIndex,
                                   moris::moris_index          aElementOwner,
                                   mtk::Cell_Info const *      aCellInfo,
                                   moris::Cell< mtk::Vertex* > aVertices):
                           mCellInfo(aCellInfo),
                           mCellId(aElementId),
                           mCellInd(aElementIndex),
                           mCellOwner(aElementOwner),
                           mCellVertices(aVertices)
                           {}

    // ----------------------------------------------------------------------------------
    // Cell get functions
    // ----------------------------------------------------------------------------------
    Matrix< DDRMat >
    Cell_XTK_No_CM::get_vertex_coords() const
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

    moris::Cell<moris::mtk::Vertex const *>
    Cell_XTK_No_CM::get_vertices_on_side_ordinal(moris::moris_index aSideOrdinal) const
    {
        moris::Cell< moris::mtk::Vertex* > tVertices = this->get_vertex_pointers();

        moris::Matrix<moris::IndexMat> tNodeOrdsOnSide = mCellInfo->get_node_to_facet_map(aSideOrdinal);

        moris::Cell<moris::mtk::Vertex const *> tVerticesOnSide(tNodeOrdsOnSide.numel());
        for(moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++)
        {
            tVerticesOnSide(i) = tVertices(tNodeOrdsOnSide(i));
        }
        return tVerticesOnSide;
    }

}



