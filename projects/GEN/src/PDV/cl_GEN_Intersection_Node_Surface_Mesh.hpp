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

namespace moris::gen
{
    class Surface_Mesh_Geometry;

    class Intersection_Node_Surface_Mesh : public Intersection_Node
    {
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
                uint                     aNodeIndex,
                const Cell< Background_Node* >& aBackgroundNodes,
                const Parent_Node&       aFirstParentNode,
                const Parent_Node&       aSecondParentNode,
                mtk::Geometry_Type       aBackgroundGeometryType,
                mtk::Interpolation_Order aBackgroundInterpolationOrder,
                Surface_Mesh_Geometry&   aInterfaceGeometry );

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

      private:

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
}
