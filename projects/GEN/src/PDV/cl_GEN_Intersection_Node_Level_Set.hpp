/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Level_Set.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node.hpp"

namespace moris::ge
{
    class Level_Set_Geometry;

    class Intersection_Node_Level_Set : public Intersection_Node
    {
      protected:
        std::weak_ptr< Level_Set_Geometry > mInterfaceGeometry;

      public:
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
                uint                                  aNodeIndex,
                const Cell< Node* >&                  aBaseNodes,
                const Parent_Node&                    aFirstParentNode,
                const Parent_Node&                    aSecondParentNode,
                real                                  aLocalCoordinate,
                mtk::Geometry_Type                    aBaseGeometryType,
                std::shared_ptr< Level_Set_Geometry > aInterfaceGeometry );

        /**
         * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
         * ancestor nodes.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
         */
        void append_dcoordinate_dadv(
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor ) override;

        /**
         * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
         *
         * @return ADV IDs
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids() override;

      private:

        /**
         * Determines if the parent nodes are intersected.
         * Used by initialize() to set mIsIntersected. Also sets mFirstParentOnInterface and mSecondParentOnInterface
         * Implementation provided here for parent class.
         *
         * @return if the parent nodes are intersected
         * @return false if there is no intersection detected
         */
        bool determine_is_intersected() override;
    };

}
