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
    class Geometry;

    class Intersection_Node_Surface_Mesh : public Intersection_Node
    {

        // variables

      public:

      protected:
        std::weak_ptr< Surface_Mesh_Geometry > mInterfaceGeometry;


      private:
        // functions

      public:
        /**
         * Constructor
         *
         * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
         * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
         * @param aFirstParentNodeIndex Index of the first parent of this node
         * @param aSecondParentNodeIndex Index of the second parent of this node
         * @param aFirstParentNodeLocalCoordinates Global coordinates of the first parent node
         * @param aSecondParentNodeLocalCoordinates Global coordinates of the second parent node
         * @param aAncestorNodeIndices Ancestor node indices
         * @param aAncestorNodeCoordinates Ancestor node global coordinates
         */
        Intersection_Node_Surface_Mesh(
                std::shared_ptr< Intersection_Node > aFirstParentNode,
                std::shared_ptr< Intersection_Node > aSecondParentNode,
                uint                                 aFirstParentNodeIndex,
                uint                                 aSecondParentNodeIndex,
                const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                const Matrix< DDUMat >&              aAncestorNodeIndices,
                const Cell< Matrix< DDRMat > >&      aAncestorNodeCoordinates );

      protected:

      private:
        /**
         * Computes the global coordinates of the intersection and the parents.
         * Used by setup() to set global coordinate member data. Implementation provided here for parent class.
         *
         * @return Matrix< DDRMat > Global location of the intersection node and its parents
         */
        Matrix< DDRMat > compute_global_coordinates() override;

        /**
         * Determines if the parent nodes are intersected.
         * Used by initialize() to set mIsIntersected. Also sets mFirstParentOnInterface and mSecondParentOnInterface
         * Implementation provided here for parent class.
         *
         * @return if the parent nodes are intersected
         * @return false if there is no intersection detected
         */
        bool determine_is_intersected(
                const Element_Intersection_Type aAncestorBasisFunction,
                const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates ) override;

        /**
         * Interpolate and return the local coordinates of this intersection node. Used to clean up constructor.
         *
         * @param aFirstNodeCoordinates Global coordinates of the first parent of this node
         * @param aSecondNodeCoordinates Global coordinates of the second parent of this node
         * @param aAncestorNodeCoordinates Global coordinates of surface mesh nodes
         * @return Local coordinate of the intersection with respect to the first parent node
         */
        real compute_local_coordinate(
                const Matrix< DDRMat >&         aFirstParentNodeCoordinates,
                const Matrix< DDRMat >&         aSecondParentNodeCoordinates,
                const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates );

        /**
         * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
         * ancestor nodes.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
         */
        void get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) override;

        /**
         * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
         *
         * @return ADV IDs
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids() override;

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
    };
}    // namespace moris::ge

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_SURFACE_MESH