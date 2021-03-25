#include "cl_MTK_Periodic_Boundary_Condition_Helper.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Cluster.hpp"
#include "cl_MTK_Double_Side_Set.hpp"

#include "fn_Parsing_Tools.hpp"
namespace moris
{
    namespace mtk
    {

        Periodic_Boundary_Condition_Helper::Periodic_Boundary_Condition_Helper( 
                                                std::shared_ptr<Mesh_Manager> aMeshManager,
                                                moris_index                   aMeshIndex,
                                                moris::ParameterList &        aParameterList):
                mMeshManager(aMeshManager),
                mMeshIndex(aMeshIndex)
        {
            std::string tMeshSideSetNames = aParameterList.get<std::string>("periodic_side_set_pair");

            string_to_cell_of_cell(aParameterList.get<std::string>("periodic_side_set_pair"), mMeshSideSetPairs);

            for(moris::uint i = 0; i < mMeshSideSetPairs.size(); i++)
            {
                for(moris::uint j = 0; j < mMeshSideSetPairs(i).size(); j++)
                {
                    std::cout<<" i = "<<i <<" | j = "<<j<<" "<<mMeshSideSetPairs(i)(j)<<std::endl;
                }
            }
        }

        void
        Periodic_Boundary_Condition_Helper::setup_periodic_boundary_conditions()
        {

            // access the integration mesh
            moris::mtk::Integration_Mesh* tIntegrationMesh = mMeshManager->get_integration_mesh(mMeshIndex);

            // get the first side set from the mesh
            moris::mtk::Set * tSet1 = tIntegrationMesh->get_set_by_name(mMeshSideSetPairs(0)(0));

            // get clusters in the set
            moris::Cell<moris::mtk::Cluster const *> tSetClusters = tSet1->get_clusters_on_set();

            // get the integration cells in the cluster
            moris::Cell<moris::mtk::Cell const *> const & tCellsInCluster = tSetClusters(0)->get_primary_cells_in_cluster( );

            std::cout<<"Integration Cell Id Set 1 = "<<tCellsInCluster(0)->get_id()<<std::endl;
            moris::print(tCellsInCluster(0)->get_vertex_ids(),"Vertex Ids");
            moris::print(tCellsInCluster(0)->get_vertex_coords(),"Vertex Coordinates");
    
            moris::Matrix<moris::IdMat> tVertexIdsInCluster = tSetClusters(0)->get_vertex_ids_in_cluster();
            moris::print(tVertexIdsInCluster,"tVertexIdsInCluster");


            // get the second side set from the mesh
            moris::mtk::Set * tSet2 = tIntegrationMesh->get_set_by_name(mMeshSideSetPairs(0)(1));

            // get clusters in the set
            moris::Cell<moris::mtk::Cluster const *> tSetClusters2 = tSet2->get_clusters_on_set();

            // get the integration cells in the cluster
            moris::Cell<moris::mtk::Cell const *> const & tCellsInCluster2 = tSetClusters2(0)->get_primary_cells_in_cluster( );

            std::cout<<"Integration Cell Id Set 2 = "<<tCellsInCluster2(0)->get_id()<<std::endl;
            moris::print(tCellsInCluster2(0)->get_vertex_ids(),"Vertex Ids");
            moris::print(tCellsInCluster2(0)->get_vertex_coords(),"Vertex Coordinates");

            // implement criteria to determine the periodic pairs
            // HARDCODED FOR NOW
           moris::Cell<moris::mtk::Vertex const *> tVertexPairing(4);
           tVertexPairing(0) = &tIntegrationMesh->get_mtk_vertex(tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id(3,EntityRank::NODE));  // paired with node 1
           tVertexPairing(1) = &tIntegrationMesh->get_mtk_vertex(tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id(9,EntityRank::NODE));  // paired with node 7
           tVertexPairing(2) = &tIntegrationMesh->get_mtk_vertex(tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id(12,EntityRank::NODE)); // paired with node 10
           tVertexPairing(3) = &tIntegrationMesh->get_mtk_vertex(tIntegrationMesh->get_loc_entity_ind_from_entity_glb_id(6,EntityRank::NODE));  // paired with node 4


            // construct the double side cluster
            Double_Side_Cluster* tDoubleSideCluster = new Double_Side_Cluster(tSetClusters(0),tSetClusters2(0),tVertexPairing);
            
            std::string                  tDoubleSideSetName     = "Periodic";
            moris::Cell<Cluster const *> tDoubleSideSetClusters = {tDoubleSideCluster};
            moris::Matrix<IndexMat>      tColors                = {{0}};
            uint                         tSpatialDim            = tIntegrationMesh->get_spatial_dim();

            // Construct the double side set        
            mtk::Set * tDblSideSet = new moris::mtk::Double_Side_Set(tDoubleSideSetName, tDoubleSideSetClusters, tColors, tSpatialDim);

            // add double sided periodic boundary condition to the integration mesh
            tIntegrationMesh->add_double_side_set(tDblSideSet);

        }
    }
}
