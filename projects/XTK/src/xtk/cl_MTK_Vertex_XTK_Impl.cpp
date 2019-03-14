/*
 * cl_MTK_Vertex_XTK_Impl.cpp
 *
 *  Created on: Mar 8, 2019
 *      Author: doble
 */

#include "cl_MTK_Vertex_XTK_Impl.hpp"
#include "cl_XTK_Background_Mesh.hpp"

namespace moris
{
namespace mtk
{
Matrix< DDRMat >
Vertex_XTK::get_coords() const
{
    MORIS_ASSERT(mBackgroundMeshPtr != nullptr || mBackgroundMeshVertex !=nullptr,"Background Mesh Pointer and Background Vertex pointer is null in XTK vertex");
    MORIS_ASSERT(!(mBackgroundMeshPtr == nullptr && mBackgroundMeshVertex ==nullptr),"Both pointers are not null");
    if(mBackgroundMeshPtr !=nullptr )
    {
        return mBackgroundMeshPtr->get_selected_node_coordinates_loc_inds({{mVertexIndex}});
    }
    else if(mBackgroundMeshVertex!=nullptr)
    {
        return mBackgroundMeshVertex->get_coords();
    }
    else
    {
        MORIS_ERROR(0,"Invalid get_coord implementation");
        return moris::Matrix<moris::DDRMat>(0,0);
    }
}
}
}
