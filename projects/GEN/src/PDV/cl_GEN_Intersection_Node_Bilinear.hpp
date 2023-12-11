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
         * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
         * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
         * @param aFirstParentNodeIndex Index of the first parent of this node
         * @param aSecondParentNodeIndex Index of the second parent of this node
         * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent node with respect to
         * the given ancestors
         * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent node with respect to
         * the given ancestors
         * @param aAncestorNodeIndices Ancestor node indices
         * @param aAncestorNodeCoordinates Ancestor node global coordinates
         * @param aInterpolationType Type of interpolation: bi or tri-linear
         * @param aInterfaceGeometry Geometry that intersects the parent to create this child
         */
        Intersection_Node_Bilinear(
                uint                                  aNodeIndex,
                const Cell< Node* >&                  aBaseNodes,
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                mtk::Geometry_Type                    aBaseGeometryType,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );

      private:

        /**
         * Gets the basis nodes that provided the field values for this level set intersection node to be created;
         * For multilinear intersection nodes, these are the background nodes.
         *
         * @return Basis nodes for interpolating sensitivities
         */
        const Cell< Basis_Node >& get_field_basis_nodes() override;

        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @return Local coordinate sensitivity
         */
        real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) override;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its first parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_first_parent() override;

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its second parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_second_parent() override;

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

        real compute_intersection_derivative( uint aAncestorIndex );
    };
}
