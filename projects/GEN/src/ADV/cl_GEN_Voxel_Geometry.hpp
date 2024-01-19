/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Voxel_Geometry.hpp
 *
 */

#pragma once

#include "cl_GEN_Geometry.hpp"
#include "cl_GEN_Voxel_Input.hpp"
#include "cl_Library_IO.hpp"

namespace moris::ge
{
    class Voxel_Geometry : public Geometry
    {
      private:
        std::shared_ptr< Voxel_Input > mVoxelInput;
        uint                           mIndex;

      public:
        /**
         * Constructor
         */
        Voxel_Geometry(
                std::shared_ptr< Voxel_Input > aVoxelField,
                uint                           aIndex );

        /**
         * A voxel geometry does not depend on ADVs, so this return false.
         *
         * @return false
         */
        virtual bool depends_on_advs() = 0;

        /**
         * Gets the geometric region of a node, based on this geometry. For a voxel geometry, the region is only positive
         * for a matching voxel index.
         *
         * @param aNodeIndex Node index
         * @param aNodeCoordinates Node coordinates
         * @return Geometric region enum
         */
        virtual Geometric_Region get_geometric_region(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aNodeCoordinates ) = 0;

        /**
         * Creates an intersection node based on the given information. The intersection node may or may not represent an intersection;
         * that is, its position may lie outside of the edge definition based on the given nodal coordinates. This information can be
         * requested from the created intersection node.
         *
         * @param aNodeIndex Node index of the new intersection node
         * @param aBaseNodes Base nodes of the element where the intersection lies
         * @param aFirstParentNode Node marking the starting point of the intersection edge
         * @param aSecondParentNode Node marking the ending point of the intersection edge
         * @param aBackgroundGeometryType Geometry type of the background element
         * @param aBackgroundInterpolationOrder Interpolation order of the background element
         * @return Voxel intersection node
         */
        virtual Intersection_Node* create_intersection_node(
                uint                     aNodeIndex,
                const Cell< Node* >&     aBaseNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder ) = 0;

        /**
         * A voxel geometry has no relevant MTK fields for remeshing, so this returns an empty vector.
         *
         * @return Empty vector
         */
        virtual Cell< std::shared_ptr< mtk::Field > > get_mtk_fields() = 0;
    };
}
