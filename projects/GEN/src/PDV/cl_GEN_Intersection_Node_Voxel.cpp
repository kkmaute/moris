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

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Voxel::Intersection_Node_Voxel(
            uint                     aNodeIndex,
            const Cell< Node* >&     aBackgroundNodes,
            const Parent_Node&       aFirstParentNode,
            const Parent_Node&       aSecondParentNode,
            mtk::Geometry_Type       aBackgroundGeometryType,
            mtk::Interpolation_Order aBackgroundInterpolationOrder,
            Voxel_Geometry&          aInterfaceGeometry )
            : Intersection_Node(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aInterfaceGeometry.compute_intersection_local_coordinate( aBackgroundNodes, aFirstParentNode, aSecondParentNode ),
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder )
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

}
