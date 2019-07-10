/*
 * cl_MTK_Cell_XTK_Impl.cpp
 *
 *  Created on: Feb 19, 2019
 *      Author: doble
 */

#include "cl_MTK_Cell_XTK_Impl.hpp"
#include "cl_XTK_Background_Mesh.hpp"
namespace moris
{
namespace mtk
{
    // ----------------------------------------------------------------------------------
    // Constructor/Deconstructor Source code
    // ----------------------------------------------------------------------------------
    Cell_XTK::Cell_XTK(moris::moris_id       aElementId,
                       moris::moris_index    aElementIndex,
                       moris::moris_index    aElementOwner,
                       moris::moris_index    aCMElementIndex,
                       xtk::Child_Mesh*      aChildMeshPtr,
                       xtk::Background_Mesh* aBackgroundMeshPtr):
                           mElementId(aElementId),
                           mElementIndex(aElementIndex),
                           mElementOwner(aElementOwner),
                           mCMElementIndex(aCMElementIndex),
                           mChildMeshPtr(aChildMeshPtr),
                           mBackgroundMeshPtr(aBackgroundMeshPtr)
                           {}

    moris::Cell< Vertex* >
    Cell_XTK::get_vertex_pointers() const
    {
        Matrix< IndexMat > tVertexIndices = this->get_vertex_inds();
        moris::Cell< Vertex* > tVertices(tVertexIndices.numel());

        for(moris::uint  i = 0; i < tVertices.size(); i++)
        {
            tVertices(i) = &mBackgroundMeshPtr->get_mtk_vertex(tVertexIndices(i));
        }
        return tVertices;
    }

    // ----------------------------------------------------------------------------------
    // Cell get functions
    // ----------------------------------------------------------------------------------
    Matrix< DDRMat >
    Cell_XTK::get_vertex_coords() const
    {
        return mBackgroundMeshPtr->get_selected_node_coordinates_loc_inds(this->get_vertex_inds());
    }

}
}


