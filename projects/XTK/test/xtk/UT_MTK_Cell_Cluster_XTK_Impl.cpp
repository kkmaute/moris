/*
 * cl_MTK_Cell_Cluster_XTK_Impl.cpp
 *
 *  Created on: Apr 26, 2019
 *      Author: doble
 */


#include "catch.hpp"

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "geometry/cl_Discrete_Level_Set.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_MTK_Mesh_XTK_Impl.hpp"

#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"


namespace xtk
{

// Cell cluster of an intersected background integration cell
class Cell_Cluster_Child_XTK: public moris::mtk::Cell_Cluster
{
public:

private:
    Cut_Mesh* mCutMesh;
};
}



namespace xtk
{

TEST_CASE("XTK Cell Clusters","[MTK_CLUSTER_XTK]")
{
    if(moris::par_size()== 0 || moris::par_size()== 2)
    {
    bool tOutputEnrichmentFields = true;

    // Generate mesh from string
    std::string tMeshFileName     = "generated:2x2x2";

    // Add level set field to add onto file

    // Specify field parameters
    moris::mtk::Scalar_Field_Info<DDRMat> tLSF;
    std::string tLSFName = "lsf1";
    tLSF.set_field_name(tLSFName);
    tLSF.set_field_entity_rank(moris::EntityRank::NODE);

    // Add to mesh input field container
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    add_field_for_mesh_input(&tLSF,tFieldsInfo);

    // add to mesh data input container
    moris::mtk::MtkMeshData tSuppMeshData;
    tSuppMeshData.FieldsInfo = &tFieldsInfo;

    // Create mesh with supplementary data
    moris::mtk::Mesh* tMeshData   = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, &tSuppMeshData );


    xtk::size_t tNumNodes = tMeshData->get_num_entities(moris::EntityRank::NODE);

    moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.2);

    moris_id tIndexOfNodeId1  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1,EntityRank::NODE);
    moris_id tIndexOfNodeId3  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3,EntityRank::NODE);
    moris_id tIndexOfNodeId5  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 5,EntityRank::NODE);
    moris_id tIndexOfNodeId7  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 7,EntityRank::NODE);
    moris_id tIndexOfNodeId9  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 9,EntityRank::NODE);
    moris_id tIndexOfNodeId11 = tMeshData->get_loc_entity_ind_from_entity_glb_id(11,EntityRank::NODE);
    moris_id tIndexOfNodeId13 = tMeshData->get_loc_entity_ind_from_entity_glb_id(13,EntityRank::NODE);
    moris_id tIndexOfNodeId15 = tMeshData->get_loc_entity_ind_from_entity_glb_id(15,EntityRank::NODE);
    moris_id tIndexOfNodeId17 = tMeshData->get_loc_entity_ind_from_entity_glb_id(17,EntityRank::NODE);
    moris_id tIndexOfNodeId19 = tMeshData->get_loc_entity_ind_from_entity_glb_id(19,EntityRank::NODE);
    moris_id tIndexOfNodeId21 = tMeshData->get_loc_entity_ind_from_entity_glb_id(21,EntityRank::NODE);
    moris_id tIndexOfNodeId23 = tMeshData->get_loc_entity_ind_from_entity_glb_id(23,EntityRank::NODE);
    moris_id tIndexOfNodeId25 = tMeshData->get_loc_entity_ind_from_entity_glb_id(25,EntityRank::NODE);
    moris_id tIndexOfNodeId27 = tMeshData->get_loc_entity_ind_from_entity_glb_id(27,EntityRank::NODE);


    // Bottom face
    tLevelsetVal(tIndexOfNodeId1) = 1;
    tLevelsetVal(tIndexOfNodeId3) = 1;
    tLevelsetVal(tIndexOfNodeId7) = 1;
    tLevelsetVal(tIndexOfNodeId9) = 1;

    // Top Face
    tLevelsetVal(tIndexOfNodeId19) = 1;
    tLevelsetVal(tIndexOfNodeId21) = 1;
    tLevelsetVal(tIndexOfNodeId25) = 1;
    tLevelsetVal(tIndexOfNodeId27) = 1;

    tLevelsetVal(tIndexOfNodeId5) = 1;
    tLevelsetVal(tIndexOfNodeId11) = 1;
    tLevelsetVal(tIndexOfNodeId17) = 1;
    tLevelsetVal(tIndexOfNodeId23) = 1;
    tLevelsetVal(tIndexOfNodeId15) = 1;
    tLevelsetVal(tIndexOfNodeId13) = 1;


    tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, moris::EntityRank::NODE, tLevelsetVal);
    tMeshData->mVerbose = true;
    std::string tPrefix2 = std::getenv("MORISOUTPUT");
    std::string tMeshOutputFile2 = tPrefix2 + "/enrichment_test_10_cluster_background.e";
    tMeshData->create_output_mesh(tMeshOutputFile2);

    // geometry
    Discrete_Level_Set tLevelSetMesh(tMeshData,{tLSFName});
    Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine(tLevelSetMesh,tPhaseTable);
    tGeometryEngine.mComputeDxDp = false;

    /*
     * Setup XTK Model and tell it how to cut
     */
    size_t tModelDimension = 3;
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
    tXTKModel.mSameMesh = true;

    /*
     * Decompose
     */
    tXTKModel.decompose(tDecompositionMethods);

//    tXTKModel.unzip_interface();

    // Perform the enrichment
//    tXTKModel.perform_basis_enrichment();

    moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();


    std::string tPrefix = std::getenv("MORISOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/enrichment_test_10_cluster.e";
    tCutMeshData->create_output_mesh(tMeshOutputFile);



    delete tMeshData;
    delete tCutMeshData;
    }


//    // Tetrathedral cells in material phase 1
//    CellTopology     tPhase0ChildTopo  = CellTopology::TET4;
//    Matrix<IndexMat> tCellIdsPhase0    = {{6, 8, 10, 12, 14, 16, 17, 18, 20, 31, 32, 33, 42, 43, 44, 53, 54, 55, 62, 63, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84}};
//    Matrix<IndexMat> tCellToNodePhase0 = {{29, 14, 31, 32},{13, 28, 30, 33},{35, 16, 36, 37},{14, 31, 32, 34},{16, 36, 37, 38},{39, 15, 40, 41},{13, 30, 40, 43},{15, 40, 41, 42},{30, 13, 33, 43},{13, 14, 31, 29},{28, 13, 31, 29},{13, 28, 31, 30},{16, 14, 35, 36},{14, 31, 34, 36},{14, 34, 35, 36},{15, 16, 39, 40},{16, 36, 38, 40},{16, 38, 39, 40},{15, 13, 40, 43},{15, 40, 42, 43},{13, 30, 31, 44},{14, 13, 31, 44},{25, 13, 14, 44},{14, 31, 36, 44},{16, 14, 36, 44},{25, 14, 16, 44},{16, 36, 40, 44},{15, 16, 40, 44},{25, 16, 15, 44},{30, 13, 40, 44},{13, 15, 40, 44},{13, 25, 15, 44}};
//    // Tetrathedral cells in material phase 1
//    CellTopology tPhase1ChildTopo = CellTopology::TET4;
//    Matrix<IndexMat> tCellToNodeGhost0 = {{21, 27, 31, 30},{17, 18, 21, 27},{31, 27, 34, 36},{18, 20, 22, 27},{36, 27, 38, 40},{20, 19, 23, 27},{17, 24, 19, 27},{30, 27, 31, 44},{31, 27, 36, 44},{36, 27, 40, 44},{27, 30, 40, 44},{17, 26, 18, 27},{18, 26, 20, 27},{20, 26, 19, 27},{17, 19, 26, 27},{21, 28, 30, 31},{21, 28, 31, 29},{21, 29, 31, 32},{27, 21, 31, 32},{18, 21, 27, 32},{28, 21, 30, 33},{21, 27, 30, 33},{21, 17, 27, 33},{27, 22, 34, 36},{34, 22, 35, 36},{22, 35, 36, 37},{27, 22, 36, 37},{20, 22, 27, 37},{31, 27, 32, 34},{27, 18, 32, 34},{27, 22, 18, 34},{27, 23, 38, 40},{38, 23, 39, 40},{36, 27, 37, 38},{27, 20, 37, 38},{27, 23, 20, 38},{23, 39, 40, 41},{27, 23, 40, 41},{19, 23, 27, 41},{27, 24, 42, 43},{30, 27, 40, 43},{40, 27, 42, 43},{40, 27, 41, 42},{27, 19, 41, 42},{27, 24, 19, 42},{27, 30, 33, 43},{17, 27, 33, 43},{24, 27, 17, 43}};    Matrix<IndexMat> tCellIdsGhost0 = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};
//    Matrix<IndexMat> tCellIdsGhost0    = {{5, 7, 9, 11, 13, 15, 19, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46, 47, 48, 49, 50, 51, 52, 56, 57, 58, 59, 60, 61, 64, 65, 66, 67, 68, 69, 70, 71, 72}};
//
//
//

}

}
