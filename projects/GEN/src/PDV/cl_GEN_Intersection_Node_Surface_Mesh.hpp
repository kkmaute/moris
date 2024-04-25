/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Surface_Mesh.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node.hpp"

// Forward declare SDF Facet
namespace moris::sdf
{
    class Facet;
}

namespace moris::gen
{
    class Surface_Mesh_Geometry;

    class Intersection_Node_Surface_Mesh : public Intersection_Node
    {
      private:
        sdf::Facet* mParentFacet;    // Pointer to the facet that intersected the edge to create this intersection node

      protected:
        Surface_Mesh_Geometry& mInterfaceGeometry;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBackgroundNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry (surface mesh)
         */
        Intersection_Node_Surface_Mesh(
                uint                              aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                real                              aLocalCoordinate,
                sdf::Facet*                       aParentFacet,
                mtk::Geometry_Type                aBackgroundGeometryType,
                mtk::Interpolation_Order          aBackgroundInterpolationOrder,
                Surface_Mesh_Geometry&            aInterfaceGeometry );

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

      private:
        /**
         * Recomputes the rotation matrix for this intersection node and returns it
         *
         */
        Matrix< DDRMat > get_rotation_matrix() const;

        /**
         * Computes the sensitivity of the local coordinate of a parent edge with respect to the facet vertices
         *
         * @return Vector< real > Sensitivities of the local coordinate with respect to the facet vertices. Size <Object dimension x number of vertices>
         */
        Matrix< DDRMat >
        compute_dxi_dfacet() const;

        /**
         * Gets the sensitivities of this node's global coordinates with respect to the ADVs which affect one of the
         * ancestor nodes.
         *
         * @param aCoordinateSensitivities Coordinate sensitivities matrix that gets appended to
         * @param aSensitivityFactor Matrix factor to scale this node's sensitivities based on a calling child's position and orientation.
         * This should be set to identity matrix of number of dimensions for any calls to this function outside of another intersection node.
         */
        void append_dcoordinate_dadv( Matrix< DDRMat >& aCoordinateSensitivities, const Matrix< DDRMat >& aSensitivityFactor ) const override;

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Gets the IDs of ADVs which one of the ancestors of this intersection node depends on.
         *
         * @return ADV IDs
         */
        Matrix< DDSMat > get_coordinate_determining_adv_ids() const override;

        //--------------------------------------------------------------------------------------------------------------
    };
}    // namespace moris::gen
