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

namespace moris::gen
{
    class Level_Set_Geometry;

    class Intersection_Node_Level_Set : public Intersection_Node
    {
      protected:
        Level_Set_Geometry& mInterfaceGeometry;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBackgroundElement The background element that this floating node is created in
         * @param aNodeManager Node manager needed to access information about other GEN nodes
         * @param aFirstParentNode First parent node - a GEN node that this intersection node is created between
         * @param aSecondParentNode Second parent node - a GEN node that this intersection node is created between
         * @param aInterfaceGeometry Interface geometry (level set)
         */
        Intersection_Node_Level_Set(
                uint                aNodeIndex,
                const mtk::Cell&    aBackgroundElement,
                const Node_Manager& aNodeManager,
                const Parent_Node&  aFirstParentNode,
                const Parent_Node&  aSecondParentNode,
                Level_Set_Geometry& aInterfaceGeometry );

        /**
         * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
         * nodes it depends on. Appends to a given sensitivity  matrix.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
         */
        void append_dcoordinate_dadv(
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor ) const override;

        /**
         * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
         *
         * @return ADV IDs
         */
        Vector< sint > get_coordinate_determining_adv_ids() const override;

      protected:
        /**
         * Gets the geometry that this intersection node was created on its interface.
         *
         * @return Geometry shared pointer
         */
        Geometry& get_interface_geometry() override;

        /**
         * Gets the geometry that this intersection node was created on its interface (const version)
         *
         * @return Const geometry reference
         */
        const Geometry& get_interface_geometry() const override;

        /**
         * Gets if the sensitivities for this node should be geometrically consistent with the mesh (true) or normal to the geometry isocontour (false)
         *
         * @return Sensitivity flag
         */
        bool compute_sensitivity_along_edges() const override;

      private:
        /**
         * Gets the sensitivity of this node's local coordinate within its parent edge with respect to the field
         * values on each of its ancestors.
         *
         * @param aAncestorIndex Ancestor index
         * @return Local coordinate sensitivity
         */
        virtual real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const = 0;
    };
}    // namespace moris::gen
