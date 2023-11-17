/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MTK_Cell_Cluster_XTK_Impl.cpp
 *
 */

#include "catch.hpp"

#include "cl_MTK_Cell_Cluster.hpp"
#include "cl_XTK_Cut_Mesh.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"

#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"

#include "cl_GEN_Mesh_Field.hpp"

namespace xtk
{

TEST_CASE("XTK Cell Clusters","[MTK_CLUSTER_XTK]")
{
    if(moris::par_size()== 1 || moris::par_size()== 2)
    {
        // Generate mesh from string
        std::string tMeshFileName     = "generated:2x2x2";

        // Add level set field to add onto file
        // Specify field parameters
        moris::mtk::Scalar_Field_Info<DDRMat> tLSF;
        std::string tLSFName = "lsf1";
        tLSF.set_field_name(tLSFName);
        tLSF.set_field_entity_rank(mtk::EntityRank::NODE);

        // Add to mesh input field container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        add_field_for_mesh_input(&tLSF,tFieldsInfo);

        // add to mesh data input container
        moris::mtk::MtkMeshData tSuppMeshData;
        tSuppMeshData.FieldsInfo = &tFieldsInfo;

        // Create mesh with supplementary data
        moris::mtk::Interpolation_Mesh* tMeshData   = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tSuppMeshData );

        xtk::size_t tNumNodes = tMeshData->get_num_entities(mtk::EntityRank::NODE);

        moris::Matrix<moris::DDRMat> tLevelsetVal(tNumNodes,1,-1.2);

        moris_id tIndexOfNodeId1  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId3  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId5  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 5,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId7  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 7,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId9  = tMeshData->get_loc_entity_ind_from_entity_glb_id( 9,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId11 = tMeshData->get_loc_entity_ind_from_entity_glb_id(11,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId13 = tMeshData->get_loc_entity_ind_from_entity_glb_id(13,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId15 = tMeshData->get_loc_entity_ind_from_entity_glb_id(15,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId17 = tMeshData->get_loc_entity_ind_from_entity_glb_id(17,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId19 = tMeshData->get_loc_entity_ind_from_entity_glb_id(19,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId21 = tMeshData->get_loc_entity_ind_from_entity_glb_id(21,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId23 = tMeshData->get_loc_entity_ind_from_entity_glb_id(23,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId25 = tMeshData->get_loc_entity_ind_from_entity_glb_id(25,mtk::EntityRank::NODE);
        moris_id tIndexOfNodeId27 = tMeshData->get_loc_entity_ind_from_entity_glb_id(27,mtk::EntityRank::NODE);

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
        tLevelsetVal(tIndexOfNodeId5)  = 1;
        tLevelsetVal(tIndexOfNodeId11) = 1;
        tLevelsetVal(tIndexOfNodeId17) = 1;
        tLevelsetVal(tIndexOfNodeId23) = 1;
        tLevelsetVal(tIndexOfNodeId15) = 1;
        tLevelsetVal(tIndexOfNodeId13) = 1;

        tMeshData->add_mesh_field_real_scalar_data_loc_inds(tLSFName, mtk::EntityRank::NODE, tLevelsetVal);
        tMeshData->mVerbose = false;
        std::string tMeshOutputFile2 = "./xtk_exo/xtk_cell_cluster_bm.e";
        tMeshData->create_output_mesh(tMeshOutputFile2);

        // geometry
        Cell< std::shared_ptr< ge::Level_Set_Geometry > > tGeometry( 1 );
        auto tField = std::make_shared<moris::ge::Mesh_Field >( tMeshData, tLSFName );
        tGeometry( 0 ) = std::make_shared< ge::Level_Set_Geometry >( tField );

        moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

        /*
         * Setup XTK Model and tell it how to cut
         */
        size_t tModelDimension = 3;
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,Subdivision_Method::C_HIERARCHY_TET4};
        Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
        tXTKModel.mVerbose  = false;
        /*
         * Decompose
         */
        tXTKModel.decompose(tDecompositionMethods);

        // Perform the enrichment
        tXTKModel.perform_basis_enrichment(mtk::EntityRank::NODE);

        delete tMeshData;
    }
}
}

