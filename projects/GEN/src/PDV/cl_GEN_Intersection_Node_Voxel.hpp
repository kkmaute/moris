/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Voxel.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node.hpp"

namespace moris::ge
{
    // Forward declare voxel geometry
    class Voxel_Geometry;

    class Intersection_Node_Voxel : public Intersection_Node
    {
      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor
         * @param aBaseNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry (voxel)
         */
        Intersection_Node_Voxel(
                uint                              aNodeIndex,
                const Cell< Node* >&              aBaseNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                mtk::Geometry_Type                aBackgroundGeometryType,
                mtk::Interpolation_Order          aBackgroundInterpolationOrder,
                std::shared_ptr< Voxel_Geometry > aInterfaceGeometry );

      private:

        /**
         * Gets the local coordinate of this intersection node based on the voxel geometry and
         * given parent node information.
         *
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aInterfaceGeometry Voxel geometry
         * @return Local coordinate along the parent edge
         */
        static real get_local_coordinate(
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                std::shared_ptr< Voxel_Geometry > aInterfaceGeometry );
    };
}
