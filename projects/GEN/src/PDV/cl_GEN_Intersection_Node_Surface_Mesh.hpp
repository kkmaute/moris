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
        Surface_Mesh_Geometry& mInterfaceGeometry;

      public:
        /**
         * Constructor
         *
         * @param aNodeIndex This node's index on the processor if it is admitted
         * @param aBaseNodes Background nodes of the element where this node resides
         * @param aFirstParentNode First parent node information
         * @param aSecondParentNode Second parent node information
         * @param aBackgroundGeometryType Background element geometry type
         * @param aBackgroundInterpolationOrder Background element interpolation order
         * @param aInterfaceGeometry Interface geometry (surface mesh)
         */
        Intersection_Node_Surface_Mesh(
                uint                     aNodeIndex,
                const Cell< Node* >&     aBaseNodes,
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
        virtual const Geometry& get_interface_geometry() const override;

      private:

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Computes the coordinate rotation matrix to move the x axis to point from the first parent node to the second parent node.
         * Used to raycast to determine local coordinate of intersection node.
         *
         * @return Matrix< DDRMat > direction cosine matrix used to rotate the coordinate frame for raycast
         */
        void transform_surface_mesh_to_local_coordinate(
                const Parent_Node&     aFirstParentNode,
                const Parent_Node&     aSecondParentNode,
                Surface_Mesh_Geometry& aInterfaceGeometry,
                uint&                  aRotationAxis );

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
                const Parent_Node&     aFirstParentNode,
                const Parent_Node&     aSecondParentNode,
                Surface_Mesh_Geometry& aInterfaceGeometry );

        //--------------------------------------------------------------------------------------------------------------

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
}    // namespace moris::ge

#endif    // MORIS_CL_GEN_INTERSECTION_NODE_SURFACE_MESH