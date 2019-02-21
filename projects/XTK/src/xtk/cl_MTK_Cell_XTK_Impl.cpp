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
    XTK_Cell::XTK_Cell(moris::moris_id       aElementId,
                       moris::moris_index    aElementIndex,
                       moris::moris_index    aCMElementIndex,
                       xtk::Child_Mesh*      aChildMeshPtr,
                       xtk::Background_Mesh* aBackgroundMeshPtr):
                           mElementId(aElementId),
                           mElementIndex(aElementIndex),
                           mCMElementIndex(aCMElementIndex),
                           mChildMeshPtr(aChildMeshPtr),
                           mBackgroundMeshPtr(aBackgroundMeshPtr)
                           {}
    // ----------------------------------------------------------------------------------
    // Cell get functions
    // ----------------------------------------------------------------------------------
    Matrix< DDRMat >
    XTK_Cell::get_vertex_coords() const
    {
        return mBackgroundMeshPtr->get_selected_node_coordinates_loc_inds(this->get_vertex_inds());
    }

}
}


