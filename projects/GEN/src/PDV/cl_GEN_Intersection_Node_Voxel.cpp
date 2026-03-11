/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Voxel.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Voxel.hpp"
#include "cl_GEN_Voxel_Geometry.hpp"
#include "cl_GEN_Parent_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Voxel::Intersection_Node_Voxel(
            uint                aNodeIndex,
            const mtk::Cell&    aBackgroundElement,
            const Node_Manager& aNodeManager,
            const Parent_Node&  aFirstParentNode,
            const Parent_Node&  aSecondParentNode,
            Voxel_Geometry&     aInterfaceGeometry )
            : Intersection_Node(
                      aNodeIndex,
                      aBackgroundElement,
                      aNodeManager,
                      aFirstParentNode,
                      aSecondParentNode,
                      aInterfaceGeometry.compute_intersection_local_coordinate( aBackgroundElement, aFirstParentNode, aSecondParentNode ) )
            , mInterfaceGeometry( aInterfaceGeometry )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    Geometry& Intersection_Node_Voxel::get_interface_geometry()
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Geometry& Intersection_Node_Voxel::get_interface_geometry() const
    {
        return mInterfaceGeometry;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Voxel::get_dxi_dcoordinate_first_parent() const
    {
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------


    Matrix< DDRMat >
    Intersection_Node_Voxel::get_dxi_dcoordinate_second_parent() const
    {
        return { {} };
    }

}    // namespace moris::gen
