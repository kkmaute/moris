/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_MDL_Mesh_Model_Helper.cpp
 *
 */

#include "catch.hpp"

#include "paths.hpp"

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

#include "cl_Matrix.hpp"        //LINALG
#include "linalg_typedefs.hpp"
#include "fn_equal_to.hpp" // ALG/src

#include "fn_norm.hpp"

#define protected public
#define private   public
#include "cl_MDL_Mesh_Model_Helper.hpp"
#undef protected
#undef private

namespace moris
{
    namespace mdl
    {
        TEST_CASE( "Mesh_Model_Helper_Test", "[moris],[mdl],[Mesh_Model_Helper]" )
        {
            uint p_rank = moris::par_rank();
            uint p_size = moris::par_size();

            if( p_size == 1 ) // specify it is a serial test only
            {
                std::string tPrefix = moris::get_base_moris_dir();
                std::string tMeshFileName = tPrefix + "projects/FEM/MDL/test/data/2_Blocks_1x2x1.g";
                std::cout<<"Mesh input name = "<<tMeshFileName<<std::endl;

                moris::mtk::Scalar_Field_Info<DDRMat> tNodeField;
                std::string tFieldName = "Temp_Field";
                tNodeField.set_field_name( tFieldName );
                tNodeField.set_field_entity_rank( mtk::EntityRank::NODE );

                // Initialize field information container
                moris::mtk::MtkFieldsInfo tFieldsInfo;

                // Place the node field into the field info container
                add_field_for_mesh_input( &tNodeField, tFieldsInfo );

                // Declare some supplementary fields
                mtk::MtkMeshData tMeshData;
                tMeshData.FieldsInfo = &tFieldsInfo;

                // construct the mesh data
                mtk::Interpolation_Mesh* tInterpMesh = mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, &tMeshData );
                mtk::Integration_Mesh*   tIntegMesh  = mtk::create_integration_mesh_from_interpolation_mesh( mtk::MeshType::STK, tInterpMesh );

                // place the pair in mesh manager
                mtk::Mesh_Manager tMeshManager;
                tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );

                //====================================================================================
                Mesh_Model_Helper tMeshModelHelper( &tMeshManager, 0 );

                tMeshModelHelper.mColorListBlock.resize( 2 );
                tMeshModelHelper.mColorListBlock( 0 ).set_size( 2, 1 );
                tMeshModelHelper.mColorListBlock( 0 )( 0, 0 ) = 0;
                tMeshModelHelper.mColorListBlock( 0 )( 1, 0 ) = 1;
                tMeshModelHelper.mColorListBlock( 1 ).set_size( 1, 1, 1 );

                tMeshModelHelper.mColorListSideSet.resize( 2 );

                tMeshModelHelper.compute_unique_node_lists();

                delete tInterpMesh;
                delete tIntegMesh;
            }
        }

    }/* namespace mdl */
}/* namespace moris */

