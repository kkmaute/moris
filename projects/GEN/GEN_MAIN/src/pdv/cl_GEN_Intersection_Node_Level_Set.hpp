/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Level_Set.hpp
 *
 */

#ifndef MORIS_CL_GEN_INTERSECTION_NODE_LEVEL_SET
#define MORIS_CL_GEN_INTERSECTION_NODE_LEVEL_SET

#include "cl_GEN_Intersection_Node.hpp"

namespace moris
{
    namespace ge
    {
        class Geometry;

        //------------------------------------------------------------------------------

        class Intersection_Node_Level_Set : public Intersection_Node
        {
            //------------------------------------------------------------------------------

          protected:
            std::weak_ptr< Geometry > mInterfaceGeometry;

            //------------------------------------------------------------------------------

          public:
            //------------------------------------------------------------------------------

            /**
             * Constructor
             *
             * @param aLocalCoordinate Local coordinate inside of parent edge
             * @param aFirstParentNode First parent node if it is also an intersection node, otherwise nullptr
             * @param aSecondParentNode Second parent node if it is also an intersection node, otherwise nullptr
             * @param aFirstParentNodeIndex Index of the first parent of this node
             * @param aSecondParentNodeIndex Index of the second parent of this node
             * @param aFirstParentNodeLocalCoordinates Local coordinates of the first parent of this node
             * @param aSecondParentNodeLocalCoordinates Local coordinates of the second parent of this node
             * @param aAncestorNodeIndices Node indices of the ancestors of this intersection node
             * @param aAncestorNodeCoordinates Coordinates of the ancestors of this intersection node
             * @param aAncestorBasisFunction Basis function of the ancestor topology
             * @param aInterfaceGeometry Geometry that intersects the parent to create this node
             * @param aIsocontourThreshold Threshold for determining the intersection location of this node
             * @param aIsocontourTolerance Tolerance for determining interface parent nodes based on geometry value
             * @param aIntersectionTolerance Tolerance for determining interface parent nodes with intersection distance
             * @param aDetermineIsIntersected Flag whether determination of intersection should be performed
             */
            Intersection_Node_Level_Set(
                    real                                 aLocalCoordinate,
                    std::shared_ptr< Intersection_Node > aFirstParentNode,
                    std::shared_ptr< Intersection_Node > aSecondParentNode,
                    moris_index                          aFirstParentNodeIndex,
                    moris_index                          aSecondParentNodeIndex,
                    const Matrix< DDRMat >&              aFirstParentNodeLocalCoordinates,
                    const Matrix< DDRMat >&              aSecondParentNodeLocalCoordinates,
                    Matrix< DDUMat >                     aAncestorNodeIndices,
                    Cell< Matrix< DDRMat > >             aAncestorNodeCoordinates,
                    const Element_Intersection_Type      aAncestorBasisFunction,
                    std::shared_ptr< Geometry >          aInterfaceGeometry );

            //------------------------------------------------------------------------------

            /**
             * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
             * ancestor nodes.
             *
             * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
             * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
             * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
             */
            void get_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) override;

            //------------------------------------------------------------------------------

            /**
             * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
             *
             * @return ADV IDs
             */
            Matrix< DDSMat > get_coordinate_determining_adv_ids() override;

            //------------------------------------------------------------------------------

          private:
            //------------------------------------------------------------------------------

            /**
             * Computes the global coordinates of the intersection and the parents.
             * Used by setup() to set global coordinate member data. Implementation provided by child class.
             *
             * @return Matrix< DDRMat > Global location of the intersection node and its parents
             */
            virtual Matrix< DDRMat > compute_global_coordinates() = 0;

            //------------------------------------------------------------------------------

            /**
             * Compute the difference between the phi value of the first parent node and
             * the isocontour threshold of the intersecting geometry. Implementation provided by child class.
             *
             * @param aAncestorBasisFunction the basis function type of the ancestor nodes
             * @param aParentNodeLocalCoordinates the parent node whose difference from the threshold that should be comptued
             * @param aParentNodeIndex the index of the parent node whose difference will be computed
             *
             * @return Phi value of first parent minus the level set value that determines the interface
             */
            virtual real compute_diff_from_threshold(
                    const Element_Intersection_Type aAncestorBasisFunction,
                    const Matrix< DDRMat >&         aParentNodeLocalCoordinates,
                    moris_index                     aParentNodeIndex ) = 0;

            //------------------------------------------------------------------------------

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

            //------------------------------------------------------------------------------

            /**
             * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
             * values on each of its ancestors.
             *
             * @param aAncestorIndex Ancestor index
             * @return Local coordinate sensitivity
             */
            virtual real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) = 0;

            //------------------------------------------------------------------------------

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its first parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix< DDRMat > get_dxi_dcoordinate_first_parent() = 0;

            //------------------------------------------------------------------------------

            /**
             * Gets the sensitivities of this node's local coordinate within its parent edge with respect to the global
             * coordinate values of its second parent.
             *
             * @return Local coordinate sensitivity
             */
            virtual Matrix< DDRMat > get_dxi_dcoordinate_second_parent() = 0;

            //------------------------------------------------------------------------------

        };    // class ge::Intersection_Node_Level_Set

        //------------------------------------------------------------------------------

    }    // namespace ge

}    // namespace moris

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_LEVEL_SET