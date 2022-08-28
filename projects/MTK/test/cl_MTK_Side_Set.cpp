/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Side_Set.cpp
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
#include "cl_MTK_Mesh.hpp"
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

namespace moris
{
namespace mtk
{
TEST_CASE("MTK Side","[MTK],[MTK_Side]")
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

        REQUIRE(tIntegMesh->get_num_side_sets() == 14);

        mtk::Set * tSideSet1 = tIntegMesh->get_set_by_index( 2 );
        mtk::Set * tSideSet2 = tIntegMesh->get_set_by_index( 3 );

        Matrix< IndexMat > tSideOrdinal1= tSideSet1->get_clusters_by_index( 0 )
                                                   ->get_cell_side_ordinals();

        Matrix< IndexMat > tVertex1= tSideSet1->get_clusters_by_index( 0 )
                                              ->get_primary_cells_in_cluster()(0)
                                              ->get_vertices_ind_on_side_ordinal(tSideOrdinal1(0,0));

        Matrix< IndexMat > tSideOrdinal2= tSideSet2->get_clusters_by_index( 1 )
                                                   ->get_cell_side_ordinals();

        Matrix< IndexMat > tVertex2= tSideSet2->get_clusters_by_index( 1 )
                                              ->get_primary_cells_in_cluster()(0)
                                              ->get_vertices_ind_on_side_ordinal(tSideOrdinal1(0,0));

//        print(tVertex1,"tVertex1");
//        print(tVertex2,"tVertex2");

        REQUIRE(tVertex1(0,0) == 10);                   REQUIRE(tVertex2(0,0) == 0);
        REQUIRE(tVertex1(0,1) == 11);                   REQUIRE(tVertex2(0,1) == 8);
        REQUIRE(tVertex1(0,2) == 8);                    REQUIRE(tVertex2(0,2) == 9);
        REQUIRE(tVertex1(0,3) == 9);                    REQUIRE(tVertex2(0,3) == 1);

        //delete tInterpMesh;
        delete tIntegMesh;
    }
}

}
}

