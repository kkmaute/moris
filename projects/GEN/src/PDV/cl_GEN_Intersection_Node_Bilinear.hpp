/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Bilinear.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node_Level_Set.hpp"

namespace moris::ge
{
    class Intersection_Node_Bilinear : public Intersection_Node_Level_Set
    {
      private:
        Matrix< DDRMat > mFirstParentNodeParametricCoordinates;
        Matrix< DDRMat > mSecondParentNodeParametricCoordinates;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBaseNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry (level set)
         */
        Intersection_Node_Bilinear(
                uint                                  aNodeIndex,
                const Cell< Node* >&                  aBaseNodes,
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                mtk::Geometry_Type                    aBackgroundGeometryType,
                mtk::Interpolation_Order              aBackgroundInterpolationOrder,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );

      private:

        /**
         * Gets the basis nodes that provided the field values for this level set intersection node to be created;
         * For multilinear intersection nodes, these are the background nodes.
         *
         * @return Basis nodes for interpolating sensitivities
         */
        const Cell< Basis_Node >& get_field_basis_nodes() const override;

        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @return Local coordinate sensitivity
         */
        real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const override;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its first parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_first_parent() const override;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its second parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_second_parent() const override;

        /**
         * Interpolate and return the local coordinates of this intersection node. Used to clean up constructor.
         *
         * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent node with respect to
         * the given ancestors
         * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent node with respect to
         * the given ancestors
         * @param aAncestorNodeIndices Ancestor node indices
         * @param aAncestorNodeCoordinates Ancestor node coordinates
         * @param aInterfaceGeometry Geometry that intersects the parent to create this child
         * @return Local coordinates
         */
        static real compute_local_coordinate(
                const Cell< Node* >&                  aBaseNodes,
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );
    };
}