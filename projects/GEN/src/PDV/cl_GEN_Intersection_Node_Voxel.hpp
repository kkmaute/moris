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

namespace moris::gen
{
    // Forward declare voxel geometry
    class Voxel_Geometry;

    class Intersection_Node_Voxel : public Intersection_Node
    {
      private:
        Voxel_Geometry& mInterfaceGeometry;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBackgroundNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry (voxel)
         */
        Intersection_Node_Voxel(
                uint                     aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder,
                Voxel_Geometry&          aInterfaceGeometry );

      protected:

        /**
         * Gets the geometry that this intersection node was created on its interface.
         *
         * @return Geometry shared pointer
         */
        Geometry& get_interface_geometry() override;

        /**
         * Gets the geometry that this intersection node was created on its interface (const version)
         *
         * @return Const geometry reference
         */
        const Geometry& get_interface_geometry() const override;
    };
}
