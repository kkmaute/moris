/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Algorithm.hpp
 *
 */

#ifndef MORIS_CL_XTK_DECOMPOSITION_ALGORITHM_HPP_
#define MORIS_CL_XTK_DECOMPOSITION_ALGORITHM_HPP_

#include "cl_XTK_Integration_Mesh_Generator.hpp"

namespace moris::xtk
{

    enum class Decomposition_Algorithm_Type
    {
        REGULAR_TEMPLATE_NONCONFORMING,
        NODE_HEIRARCHY,
        OCTREE,
        ELEVATE_ORDER,
        DELAUNAY,
        MAX_ENUM
    };

    class Cut_Integration_Mesh;
    class Integration_Mesh_Generator;

    struct Decomposition_Data;
    struct Integration_Mesh_Generation_Data;

    class Decomposition_Algorithm
    {
      public:
        // number of cells less the number you are replace
        moris_index                                        mNumNewCells = 0;
        Vector< Vector< moris::moris_index > >             mNewCellToVertexConnectivity;    // over allocated
        Vector< moris::moris_index >                       mNewCellChildMeshIndex;          // over allocated
        Vector< moris::moris_index >                       mNewCellCellIndexToReplace;      // over allocated
        Vector< std::shared_ptr< moris::mtk::Cell_Info > > mNewCellCellInfo;                // over allocated

      public:
        Decomposition_Algorithm() {}
        virtual ~Decomposition_Algorithm() {}

        /**
         * Checks if the given cell should be operated on by this decomposition algorithm
         *
         * @param aCell The cell to check
         * @return True if the cell will be decomposed by this algorithm, false otherwise
         */
        virtual bool is_eligible( std::pair< mtk::Cell*, Vector< Decomposition_Algorithm_Type > >& aElementContext,
                Cut_Integration_Mesh*                                                              aCutIntegrationMesh,
                Integration_Mesh_Generator*                                                        aMeshGenerator ) const = 0;

        /**
         * Gets the indices of the cells in the background mesh that were decomposed by this algorithm
         */
        virtual Vector< moris_index > get_decomposed_cell_indices() = 0;

        // set of
        virtual void perform(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator );

        virtual enum Decomposition_Algorithm_Type get_algorithm_type() const = 0;

        virtual moris_index get_signature() const = 0;

        virtual bool has_geometric_independent_vertices() const = 0;

        virtual void
        perform_impl_vertex_requests(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) = 0;

        virtual void
        perform_impl_generate_mesh(
                Integration_Mesh_Generation_Data* aMeshGenerationData,
                Decomposition_Data*               aDecompositionData,
                Cut_Integration_Mesh*             aCutIntegrationMesh,
                moris::mtk::Mesh*                 aBackgroundMesh,
                Integration_Mesh_Generator*       aMeshGenerator ) = 0;
    };

}    // namespace moris::xtk

#endif
