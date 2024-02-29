/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Decomposition_Algorithm.cpp
 *
 */

#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"
#include "cl_Tracer.hpp"

namespace moris::xtk
{

    void
    Decomposition_Algorithm::perform(
            Integration_Mesh_Generation_Data* aMeshGenerationData,
            Decomposition_Data*               aDecompositionData,
            Cut_Integration_Mesh*             aCutIntegrationMesh,
            moris::mtk::Mesh*                 aBackgroundMesh,
            Integration_Mesh_Generator*       aMeshGenerator )
    {
        Tracer tTracer( "XTK", "Decomposition_Algorithm", "Generic", aMeshGenerator->verbosity_level(), 0 );

        this->perform_impl_vertex_requests( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, aMeshGenerator );

        // in parallel give all these nodes ids
        aMeshGenerator->assign_node_requests_identifiers( *aDecompositionData, aCutIntegrationMesh, aBackgroundMesh );

        // commit vertices to the mesh
        aMeshGenerator->commit_new_ig_vertices_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );

        // generate_decomposed_mesh - perform the mesh generation
        this->perform_impl_generate_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, aMeshGenerator );

        // commit the cells to the mesh
        aMeshGenerator->commit_new_ig_cells_to_cut_mesh( aMeshGenerationData, aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this );
    }
}    // namespace moris::xtk
