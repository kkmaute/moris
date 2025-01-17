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
        Vector< moris_index >    mCellsWithSurfacePoints;

      public:
        Delaunay_Subdivision_Interface( Parameter_List& aParameterList );

        ~Delaunay_Subdivision_Interface() override {}

        bool is_eligible( std::pair< mtk::Cell*, Vector< Decomposition_Algorithm_Type > >& aElementContext,
                Cut_Integration_Mesh*                                                      aCutIntegrationMesh,
                Integration_Mesh_Generator*                                                aMeshGenerator ) const override;

        Vector< moris_index > get_decomposed_cell_indices() override;

        bool has_geometric_independent_vertices() const override;

        // void perform(
        //         Integration_Mesh_Generation_Data* aMeshGenerationData,
        //         Decomposition_Data*               aDecompositionData,
        //         Cut_Integration_Mesh*             aCutIntegrationMesh,
        //         mtk::Mesh*                        aBackgroundMesh,
        //         Integration_Mesh_Generator*       aMeshGenerator ) override;

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
        // BRENDAN write documentation
        Vector< Vector< moris_index > > triangulation();

        mtk::CellTopology get_ig_cell_topology() const;
    };

}    // namespace moris::xtk
#endif
