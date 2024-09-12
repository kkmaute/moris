/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Blocks.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"

// implementations to test
#include "cl_MTK_Vertex.hpp"    //MTK
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Mesh_Data_STK.hpp"
#include "cl_MTK_Mesh_Core_STK.hpp"
#include "cl_MTK_Interpolation_Mesh_STK.hpp"
#include "cl_MTK_Integration_Mesh_STK.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
#include "cl_Communication_Tools.hpp"

namespace moris::mtk
{
TEST_CASE("MTK Blocks","[MTK],[MTK_BLOCK]")
{
    if(par_size() ==1)
    {
        std::string tPrefix = moris::get_base_moris_dir();
        std::string tMeshFileName = tPrefix + "projects/FEM/MDL/test/data/2_Blocks_1x2x1.g";

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField;
        std::string tFieldName = "Temp_Field";
        tNodeField.set_field_name( tFieldName );
        tNodeField.set_field_entity_rank( EntityRank::NODE );

        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;

        // Place the node field into the field info container
        add_field_for_mesh_input( &tNodeField, tFieldsInfo );

        // Declare some supplementary fields
        mtk::MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;

        // construct the mesh data
        mtk::Interpolation_Mesh* tInterpMesh = mtk::create_interpolation_mesh( MeshType::STK, tMeshFileName, &tMeshData );
        mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );

        REQUIRE(tIntegMesh->get_num_blocks() == 2);

        mtk::Set * tBlock1 = tIntegMesh->get_set_by_index( 0 );
        mtk::Set * tBlock2 = tIntegMesh->get_set_by_index( 1 );

//        REQUIRE(tBlock1->get_list_of_block_cell_clusters()(0,0) == 0);
//        REQUIRE(tBlock2->get_list_of_block_cell_clusters()(0,0) == 1);

        Matrix< IndexMat > tVertexId1= tBlock1->get_clusters_by_index( 0 )
                                              ->get_primary_cells_in_cluster()(0)
                                              ->get_vertex_ids();

        Matrix< IndexMat > tVertexId2= tBlock2->get_clusters_by_index( 0 )
                                               ->get_primary_cells_in_cluster()(0)
                                               ->get_vertex_ids();

        REQUIRE(tVertexId1(0,0) == 4);                   REQUIRE(tVertexId2(0,0) == 5);
        REQUIRE(tVertexId1(0,1) == 7);                   REQUIRE(tVertexId2(0,1) == 11);
        REQUIRE(tVertexId1(0,2) == 8);                   REQUIRE(tVertexId2(0,2) == 12);
        REQUIRE(tVertexId1(0,3) == 3);                   REQUIRE(tVertexId2(0,3) == 8);
        REQUIRE(tVertexId1(0,4) == 1);                   REQUIRE(tVertexId2(0,4) == 2);
        REQUIRE(tVertexId1(0,5) == 6);                   REQUIRE(tVertexId2(0,5) == 9);
        REQUIRE(tVertexId1(0,6) == 5);                   REQUIRE(tVertexId2(0,6) == 10);
        REQUIRE(tVertexId1(0,7) == 2);                   REQUIRE(tVertexId2(0,7) == 3);

        REQUIRE(tBlock1->get_interpolation_cell_geometry_type()       == mtk::Geometry_Type::HEX);
        REQUIRE(tBlock1->get_integration_cell_geometry_type()         == mtk::Geometry_Type::HEX);
        REQUIRE(tBlock1->get_interpolation_cell_interpolation_order() == mtk::Interpolation_Order::LINEAR );
        REQUIRE(tBlock1->get_integration_cell_interpolation_order()   == mtk::Interpolation_Order::LINEAR);

        //delete tInterpMesh;
        delete tIntegMesh;
    }
}

}
