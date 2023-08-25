/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Bilinear.hpp
 *
 */

#ifndef MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP

#include "cl_GEN_Intersection_Node_Level_Set.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        class Intersection_Node_Bilinear : public Intersection_Node_Level_Set
        {
          private:
            Matrix< DDRMat > mParentLocalCoordinates;

            Element_Intersection_Type mInterpolationType;

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
             * @param aAncestorNodeCoordinates Ancestor node coordinates
             * @param aInterpolationType Type of interpolation: bi or tri-linear
             * @param aInterfaceGeometry Geometry that intersects the parent to create this child
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining interface parent nodes based on geometry value
             * @param aIntersectionTolerance Tolerance for determining interface parent nodes with intersection distance
             */
            Intersection_Node_Bilinear(
                    std::shared_ptr< Intersection_Node > aFirstParentNode,
                    std::shared_ptr< Intersection_Node > aSecondParentNode,
                    uint                                 aFirstParentNodeIndex,
                    uint                                 aSecondParentNodeIndex,
                    const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                    const Matrix< DDUMat >&              aAncestorNodeIndices,
                    const Cell< Matrix< DDRMat > >&      aAncestorNodeCoordinates,
                    const Element_Intersection_Type      aInterpolationType,
                    std::shared_ptr< Geometry >          aInterfaceGeometry );

            virtual ~Intersection_Node_Bilinear();

          private:
            /**
             * Computes the global coordinates of the intersection and the parents.
             * Used by setup() to set global coordinate member data. Implementation provided for parent class.
             *
             * @return Matrix< DDRMat > Global location of the intersection node and its parents
             */
            virtual Matrix< DDRMat > compute_global_coordinates() override;

            /**
             * Computes the vector from the first parent to the second parent
             * Used by setup to set mParentVector member data. Implementation provided by child class.
             *
             * @return Matrix< DDRMat >
             */
            Matrix< DDRMat > compute_parent_vector() override;

            /**
             * Determines if the first parent is on an interface. Used by setup() to assign member data.
             * Implemented provided for parent class.
             *
             * @return true if the first parent is on the interface
             * @return false if the first parent is not on the interface
             */
            bool determine_first_parent_on_interface() override;

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @return Local coordinate sensitivity
             */
            real get_dxi_dfield_from_ancestor( uint aAncestorIndex );

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix< DDRMat > get_dxi_dcoordinate_first_parent();

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix< DDRMat > get_dxi_dcoordinate_second_parent();

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
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @return Local coordinates
             */

            real compute_local_coordinate(
                    const Matrix< DDRMat >&         aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&         aSecondParentNodeLocalCoordinates,
                    const Matrix< DDUMat >&         aAncestorNodeIndices,
                    const Cell< Matrix< DDRMat > >& aAncestorNodeCoordinates,
                    std::shared_ptr< Geometry >     aInterfaceGeometry );

            real compute_intersection_derivative( uint aAncestorIndex );
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_BILINEAR_HPP
