/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Delaunay_Subdivision_Interface.hpp
 *
 */

#ifndef MORIS_CL_XTK_Delaunay_Subdivision_Interface_HPP_
#define MORIS_CL_XTK_Delaunay_Subdivision_Interface_HPP_

#include "cl_XTK_Decomposition_Algorithm.hpp"

namespace moris::mtk
{
    class Mesh;
}

namespace moris::xtk
{
    struct Delaunay_Template
    {
        virtual Matrix< DDRMat > get_node_param_coords( uint aPoint ) const = 0;
        virtual uint             get_num_corner_nodes() const               = 0;
        virtual uint             get_num_new_template_nodes() const         = 0;
        virtual Vector< uint >   get_new_node_face_ordinal() const          = 0;
    };

    struct Delaunay_Template_QUAD4 : public Delaunay_Template
    {
      public:
        Matrix< DDRMat > get_node_param_coords( uint aPoint ) const override
        {
            switch ( aPoint )
            {
                case 0:
                    return { { -1.0, -1.0 } };
                case 1:
                    return { { 1.0, -1.0 } };
                case 2:
                    return { { 1.0, 1.0 } };
                case 3:
                    return { { -1.0, 1.0 } };
                default:
                    MORIS_ERROR( false, "Invalid point index for delaunay template." );
                    return { { 0.0, 0.0 } };
            }
        }

        uint get_num_corner_nodes() const override
        {
            return 4;
        }

        uint get_num_new_template_nodes() const override
        {
            return 0;    // No new derived nodes for QUAD4
        }

        Vector< uint > get_new_node_face_ordinal() const override
        {
            return { MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX };    // No new derived nodes for QUAD4
        }
    };

    struct Delaunay_Template_HEX8 : public Delaunay_Template
    {
      public:
        Matrix< DDRMat > get_node_param_coords( uint aPoint ) const override
        {
            switch ( aPoint )
            {
                case 0:
                    return { { -1.0, -1.0, -1.0 } };
                case 1:
                    return { { 1.0, -1.0, -1.0 } };
                case 2:
                    return { { 1.0, 1.0, -1.0 } };
                case 3:
                    return { { -1.0, 1.0, -1.0 } };
                case 4:
                    return { { -1.0, -1.0, 1.0 } };
                case 5:
                    return { { 1.0, -1.0, 1.0 } };
                case 6:
                    return { { 1.0, 1.0, 1.0 } };
                case 7:
                    return { { -1.0, 1.0, 1.0 } };
                case 8:
                    return { { 0.0, -1.0, 0.0 } };
                case 9:
                    return { { 1.0, 0.0, 0.0 } };
                case 10:
                    return { { 0.0, 1.0, 0.0 } };
                case 11:
                    return { { -1.0, 0.0, 0.0 } };
                case 12:
                    return { { 0.0, 0.0, -1.0 } };
                case 13:
                    return { { 0.0, 0.0, 1.0 } };
                default:
                    MORIS_ERROR( false, "Invalid point index for delaunay template." );
                    return { { 0.0, 0.0, 0.0 } };
            }
        }

        uint get_num_corner_nodes() const override
        {
            return 8;
        }

        uint get_num_new_template_nodes() const override
        {
            return 6;    // 1 new node per face (+ surface point nodes later on)
        }

        Vector< uint > get_new_node_face_ordinal() const override
        {
            return { MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, MORIS_UINT_MAX, 0, 1, 2, 3, 4, 5 };
        }
    };

    // -------------------------------------------------------------------------

    class Delaunay_Subdivision_Interface : public Decomposition_Algorithm
    {
      private:
        moris::gen::Geometry_Engine*         mGeometryEngine;
        Integration_Mesh_Generation_Data*    mMeshGenerationData;
        Decomposition_Data*                  mDecompositionData;
        Cut_Integration_Mesh*                mCutIntegrationMesh;
        moris::mtk::Mesh*                    mBackgroundMesh;
        Integration_Mesh_Generator*          mGenerator;
        moris::uint                          mNumTotalCells = 0;
        std::shared_ptr< Delaunay_Template > mDelaunayTemplate;

        Vector< Vector< real > > mAllSurfacePoints;

      public:
        Delaunay_Subdivision_Interface( Parameter_List& aParameterList, mtk::CellTopology aCellTopology );

        ~Delaunay_Subdivision_Interface() override {}

        Vector< moris_index > get_decomposed_cell_indices() override;

        bool has_geometric_dependent_vertices() const override;

        void
        perform_impl_vertex_requests(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) override;

        void
        perform_impl_generate_mesh(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) override;

        enum Decomposition_Algorithm_Type
        get_algorithm_type() const override;

        moris_index get_signature() const override;

        // template functions

      private:
        /**
         * Performs a constrained Delaunay triangulation using geompack3d for all child meshes with surface points.
         * Uses mAllSurfacePoints to get the points for triangulation, which must be computed first.
         *
         * @return Vector< Vector< moris_index > > - the connectivity of the new cells. Outer index is the child mesh,
         * inner index is the local triangulation indices, starting from 1.
         */
        Vector< Vector< moris_index > > triangulation();

        /**
         * Checks that the volume of the tetrahedra is positive, inplace reorders vertices if this is not the case
         *
         * @param aCoordinates Coordinates of all the vertices in a given child mesh. Flattened vector where each group of 3 is a vertex
         * @param aConnectivity Tet connectivity of all the tets in the child mesh. Flattened vector where each group of 4 is a tet. This is modified in place
         * @param aNumTetrahedra Number of tetrahedra in the child mesh, since aConnectivity is oversized
         */
        void reorder_tet_connectivity( const Vector< real >& aCoordinates, Vector< uint >& aConnectivity, uint aNumTetrahedra );

        /**
         * Computes the sign of the volume of a tetrahedron given the coordinates of its vertices.
         * NOTE: the true volume is 1/6 of the returned value
         *
         * @param aTetCoords Coordinates of the vertices of the tetrahedron. <3, 4>
         */
        static real compute_tet_volume_sign( const Matrix< DDRMat >& aTetCoords );

        mtk::CellTopology get_ig_cell_topology() const;

        /**
         * Computes the global coordinate matrix from the parametric coordinates
         *
         * @param aParametricCoordinates Parametric coordinates of all the surface points within the cell. <nPoints, nDim>
         * @param aCell The background cell that the surface points lie in
         * @return Matrix< DDRMat > Global coordinates of the surface points. <nPoints, nDim>
         */
        static Matrix< DDRMat > get_surface_point_global_coordinates( const Matrix< DDRMat >& aParametricCoordinates, const mtk::Cell& aCell );

        /**
         * Gets the number of corner nodes for a given cell type
         *
         * @param aCellType The type of cell (tri, quad, tet, hex, etc.)
         * @return the number of corner nodes
         */
        uint get_num_geometric_nodes( const mtk::Geometry_Type aCellType ) const;
    };

}    // namespace moris::xtk
#endif
