/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.hpp
 *
 */

#ifndef MORIS_CL_GEN_INTERSECTION_NODE_SURFACE_MESH
#define MORIS_CL_GEN_INTERSECTION_NODE_SURFACE_MESH

#include "cl_GEN_Intersection_Node.hpp"


namespace moris::ge
{
    class Surface_Mesh_Geometry;

    class Intersection_Node_Surface_Mesh : public Intersection_Node
    {
      protected:
        std::weak_ptr< Surface_Mesh_Geometry > mInterfaceGeometry;

      private:
        // functions

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex the index assigned to this node
         * @param aBaseNodes the background nodes that correspond to this intersection node
         * @param aFirstParentNode the first node that this intersection node lies between. can be a base node or a derived node such as another intersection node
         * @param aSecondParentNode the other node that this intersection node lies between. can be a base node or a derived node such as another intersection node
         * @param aBaseGeometryType type of collection of base nodes. QUAD for 2D and HEX for 3D
         * @param aInterfaceGeometry geometry that intersects the parents to create this intersection node
         */
        Intersection_Node_Surface_Mesh(
                uint                                     aNodeIndex,
                const Cell< Node* >&                     aBaseNodes,
                const Parent_Node&                       aFirstParentNode,
                const Parent_Node&                       aSecondParentNode,
                mtk::Geometry_Type                       aBaseGeometryType,
                std::shared_ptr< Surface_Mesh_Geometry > aInterfaceGeometry );

      protected:

      private:
        //--------------------------------------------------------------------------------------------------------------

        /**
         * Determines if the parent nodes are intersected.
         * Used by initialize() to set mIsIntersected. Also sets mFirstParentOnInterface and mSecondParentOnInterface
         * Implementation provided here for parent class.
         *
         * @return if the parent nodes are intersected
         * @return false if there is no intersection detected
         */
        bool determine_is_intersected() override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Computes the coordinate rotation matrix to move the x axis to point from the first parent node to the second parent node.
         * Used to raycast to determine local coordinate of intersection node.
         *
         * @return Matrix< DDRMat > direction cosine matrix used to rotate the coordinate frame for raycast
         */
        Matrix< DDRMat > compute_raycast_rotation();

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Interpolate and return the local coordinates of this intersection node by raycasting in a transformed coordinate axis z direction.
         * Used to clean up constructor.
         *
         * @param aFirstNodeCoordinates Global coordinates of the first parent of this node
         * @param aSecondNodeCoordinates Global coordinates of the second parent of this node
         * @param aAncestorNodeCoordinates Global coordinates of surface mesh nodes
         * @return Local coordinate of the intersection with respect to the first parent node
         */
        real compute_local_coordinate(
                const Parent_Node& aFirstParentNode,
                const Parent_Node& aSecondParentNode );

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
         * ancestor nodes.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
         */
        void append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
         *
         * @return ADV IDs
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids() override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @param aAncestorIndex Ancestor index
         * @return Local coordinate sensitivity
         */
        real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its first parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_first_parent() override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
         * coordinate values of its second parent.
         *
         * @return Local coordinate sensitivity
         */
        Matrix< DDRMat > get_dxi_dcoordinate_second_parent() override;

        //--------------------------------------------------------------------------------------------------------------
    };
}    // namespace moris::ge

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_SURFACE_MESH