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
    // -------------------------------------------------------------------------

    class Delaunay_Subdivision_Interface : public Decomposition_Algorithm
    {
      private:
        moris::gen::Geometry_Engine*      mGeometryEngine;
        Integration_Mesh_Generation_Data* mMeshGenerationData;
        Decomposition_Data*               mDecompositionData;
        Cut_Integration_Mesh*             mCutIntegrationMesh;
        moris::mtk::Mesh*                 mBackgroundMesh;
        Integration_Mesh_Generator*       mGenerator;
        moris::uint                       mNumTotalCells = 0;

        Vector< Vector< real > > mAllSurfacePoints;

      public:
        Delaunay_Subdivision_Interface( Parameter_List& aParameterList );

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

        mtk::CellTopology get_ig_cell_topology() const;

        /**
         * Computes the global coordinate matrix from the parametric coordinates
         * 
         * @param aParametricCoordinates Parametric coordinates of all the surface points within the cell. <nPoints, nDim>
         * @param aCell The background cell that the surface points lie in
         * @return Matrix< DDRMat > Global coordinates of the surface points. <nPoints, nDim>
         */
        Matrix< DDRMat > get_surface_point_global_coordinates( const Matrix< DDRMat >& aParametricCoordinates, const mtk::Cell& aCell );

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
