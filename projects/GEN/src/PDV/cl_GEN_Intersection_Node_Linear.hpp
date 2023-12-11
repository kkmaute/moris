/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Linear.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node_Level_Set.hpp"

namespace moris::ge
{
    class Field;

    class Intersection_Node_Linear : public Intersection_Node_Level_Set
    {
      public:
        /**
         * Constructor
         *
         * @param aFirstNode First parent node if it is also an intersection node, otherwise nullptr
         * @param aSecondNode Second parent node if it is also an intersection node, otherwise nullptr
         * @param aFirstNodeIndex Index of the first parent of this node
         * @param aSecondNodeIndex Index of the second parent of this node
         * @param aFirstNodeCoordinates Coordinates of the first parent of this node
         * @param aSecondNodeCoordinates Coordinates of the second parent of this node
         * @param aInterfaceGeometry Geometry that intersects the parent to create this child
         * @param aIntersectionTolerance Tolerance for determining interface parent nodes with intersection distance
         */
        Intersection_Node_Linear(
                uint                                  aNodeIndex,
                const Cell< Node* >&                  aBaseNodes,
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                mtk::Geometry_Type                    aBaseGeometryType,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );

        /**
         * Given a node index or coordinates, returns a vector of the field derivatives with respect to the nodal
         * coordinates.
         *
         * @param aField Field pointer, referenced during call from field class
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                Field*            aField,
                Matrix< DDRMat >& aSensitivities );

      private:

        /**
         * Gets the basis nodes that provided the field values for this level set intersection node to be created;
         * For a linear intersection node, these are the parent nodes.
         *
         * @return Basis nodes for interpolating sensitivities
         */
        const Cell< Basis_Node >& get_field_basis_nodes() override;

        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @param aAncestorIndex Ancestor index
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
         * @param aFirstNodeIndex Index of the first parent of this node
         * @param aSecondNodeIndex Index of the second parent of this node
         * @param aFirstNodeCoordinates Coordinates of the first parent of this node
         * @param aSecondNodeCoordinates Coordinates of the second parent of this node
         * @param aInterfaceGeometry Geometry that intersects the parent to create this node
         * @return Local coordinates
         */
        real compute_local_coordinate(
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );
    };
}
