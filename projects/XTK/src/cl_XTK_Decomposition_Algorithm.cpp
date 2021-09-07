#include "cl_XTK_Decomposition_Algorithm.hpp"
#include "cl_XTK_Integration_Mesh_Generator.hpp"

namespace xtk
{

    void 
    Decomposition_Algorithm::perform(
    Integration_Mesh_Generation_Data*  aMeshGenerationData,
    Decomposition_Data*                 aDecompositionData,
    Cut_Integration_Mesh*               aCutIntegrationMesh,
    moris::mtk::Mesh*                   aBackgroundMesh,
    Integration_Mesh_Generator*         aMeshGenerator)
    {
            this->perform_impl_vertex_requests(aMeshGenerationData,aDecompositionData,aCutIntegrationMesh,aBackgroundMesh,aMeshGenerator);

            // in parallel give all these nodes ids
            aMeshGenerator->assign_node_requests_identifiers(*aDecompositionData, aCutIntegrationMesh, aBackgroundMesh, this->get_signature());

            // commit vertices to the mesh
            aMeshGenerator->commit_new_ig_vertices_to_cut_mesh(aMeshGenerationData,aDecompositionData,aCutIntegrationMesh,aBackgroundMesh,this);

            // generate_decomposed_mesh - perform the mesh generation
            this->perform_impl_generate_mesh(aMeshGenerationData,aDecompositionData,aCutIntegrationMesh,aBackgroundMesh,aMeshGenerator);

            // commit the cells to the mesh
            aMeshGenerator->commit_new_ig_cells_to_cut_mesh(aMeshGenerationData,aDecompositionData,aCutIntegrationMesh,aBackgroundMesh, this);
    }
}