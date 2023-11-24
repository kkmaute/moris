/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Linear.hpp
 *
 */

#ifndef MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP
#define MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP

#include "cl_GEN_Intersection_Node_Level_Set.hpp"

namespace moris
{
    namespace ge
    {
        class Level_Set_Geometry;

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
                    std::shared_ptr< Intersection_Node >  aFirstNode,
                    std::shared_ptr< Intersection_Node >  aSecondNode,
                    uint                                  aFirstNodeIndex,
                    uint                                  aSecondNodeIndex,
                    const Matrix< DDRMat >&               aFirstNodeCoordinates,
                    const Matrix< DDRMat >&               aSecondNodeCoordinates,
                    std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );

            //--------------------------------------------------------------------------------------------------------------

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

            //--------------------------------------------------------------------------------------------------------------

          private:
            /**
             * Computes the global coordinates of the intersection and the parents.
             * Used by setup() to set global coordinate member data. Implementation provided for parent class.
             *
             * @return Matrix< DDRMat > Global location of the intersection node and its parents
             */
            Matrix< DDRMat > compute_global_coordinates() override;

            //--------------------------------------------------------------------------------------------------------------

            /**
             * Compute the difference between the phi value of the first parent node and
             * the isocontour threshold of the intersecting geometry. Implementation provided here for parent class.
             *
             * @param aAncestorBasisFunction the basis function type of the ancestor nodes
             * @param aParentNodeLocalCoordinates the parent node whose difference from the threshold that should be comptued
             * @param aParentNodeIndex the index of the parent node whose difference will be computed
             *
             * @return Phi value of first parent minus the level set value that determines the interface
             */
            real compute_diff_from_threshold(
                    const Element_Interpolation_Type aAncestorBasisFunction,
                    const Matrix< DDRMat >&          aParentNodeLocalCoordinates,
                    moris_index                      aParentNodeIndex ) override;

            //--------------------------------------------------------------------------------------------------------------

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @param aAncestorIndex Ancestor index
             * @return Local coordinate sensitivity
             */
            real get_dxi_dfield_from_ancestor( uint aAncestorIndex );

            //--------------------------------------------------------------------------------------------------------------

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix< DDRMat > get_dxi_dcoordinate_first_parent();

            //--------------------------------------------------------------------------------------------------------------

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            Matrix< DDRMat > get_dxi_dcoordinate_second_parent();

            //--------------------------------------------------------------------------------------------------------------

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
                    uint                                  aFirstNodeIndex,
                    uint                                  aSecondNodeIndex,
                    const Matrix< DDRMat >&               aFirstNodeCoordinates,
                    const Matrix< DDRMat >&               aSecondNodeCoordinates,
                    std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );
        };
    }    // namespace ge
}    // namespace moris

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_LINEAR_HPP
