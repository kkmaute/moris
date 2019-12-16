/*
 * cl_Parameter_Sweep.cpp
 *
 *  Created on: Dec 10, 2019
 *      Author: sonne
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

// GE include -----------------------------------
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Phase_Table.hpp"

#include "cl_GEN_Field_User_Defined.hpp"

#include "cl_GEN_Geometry_Engine.hpp"

#include "cl_GEN_Sphere.hpp"

// HMR includes ---------------------------------
#include "cl_HMR.hpp"

// MTK includes ---------------------------------
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

// XTK include ----------------------------------
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
//------------------------------------------------------------------------------


namespace moris{
namespace ge{

TEST_CASE("param_test_01","[GE],[param_sweep_01]")
{
    if(par_size()<=1)
    {
        const std::string fileName = "generated:10x10x10";
        //------------------------------------------------------------------------------
        moris::Cell< std::string > tFieldNames(10);
        tFieldNames(0) = "tSphere00";   tFieldNames(5) = "tSphere05";
        tFieldNames(1) = "tSphere01";   tFieldNames(6) = "tSphere06";
        tFieldNames(2) = "tSphere02";   tFieldNames(7) = "tSphere07";
        tFieldNames(3) = "tSphere03";   tFieldNames(8) = "tSphere08";
        tFieldNames(4) = "tSphere04";   tFieldNames(9) = "tSphere09";

        // Declare scalar node fields
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField0;
        std::string tFieldName0 = tFieldNames(0);
        tNodeField0.set_field_name(tFieldName0);
        tNodeField0.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName1 = tFieldNames(1);
        tNodeField1.set_field_name(tFieldName1);
        tNodeField1.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField2;
        std::string tFieldName2 = tFieldNames(2);
        tNodeField2.set_field_name(tFieldName2);
        tNodeField2.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField3;
        std::string tFieldName3 = tFieldNames(3);
        tNodeField3.set_field_name(tFieldName3);
        tNodeField3.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField4;
        std::string tFieldName4 = tFieldNames(4);
        tNodeField4.set_field_name(tFieldName4);
        tNodeField4.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField5;
        std::string tFieldName5 = tFieldNames(5);
        tNodeField5.set_field_name(tFieldName5);
        tNodeField5.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField6;
        std::string tFieldName6 = tFieldNames(6);
        tNodeField6.set_field_name(tFieldName6);
        tNodeField6.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField7;
        std::string tFieldName7 = tFieldNames(7);
        tNodeField7.set_field_name(tFieldName7);
        tNodeField7.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField8;
        std::string tFieldName8 = tFieldNames(8);
        tNodeField8.set_field_name(tFieldName8);
        tNodeField8.set_field_entity_rank(EntityRank::NODE);

        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField9;
        std::string tFieldName9 = tFieldNames(9);
        tNodeField9.set_field_name(tFieldName9);
        tNodeField9.set_field_entity_rank(EntityRank::NODE);
        //------------------------------------------------------------------------------
        // Initialize field information container
        moris::mtk::MtkFieldsInfo tFieldsInfo;
        // Place the node field into the field info container
        add_field_for_mesh_input( &tNodeField0, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField1, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField2, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField3, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField4, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField5, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField6, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField7, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField8, tFieldsInfo );
        add_field_for_mesh_input( &tNodeField9, tFieldsInfo );
        // Declare some supplementary fields
        mtk::MtkMeshData tMeshData;
        tMeshData.FieldsInfo = &tFieldsInfo;
        // create mesh pair
        mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, fileName, &tMeshData );
        mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK, tInterpMesh );
        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
        //------------------------------------------------------------------------------
        uint tNumNodes = tInterpMesh->get_num_nodes();
        Matrix< DDRMat > tLSVals( tNumNodes,1 );
        moris::Cell<moris::mtk::Vertex const *> tAllVerices = tInterpMesh->get_all_vertices();

        Matrix< DDRMat > tAllCoords( tNumNodes, 3 );
        for(uint i=0; i<tNumNodes; i++ )
        {
            tAllCoords.set_row( i,tAllVerices(i)->get_coords() );
        }

        real tXcenter = 5.0;
        real tYcenter = 5.0;
        real tZcenter = 5.0;

        Matrix< DDRMat > tAllRadii{{ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 }};
        //------------------------------------------------------------------------------
        for(uint iIter=0; iIter<10; iIter++)
        {
            real tRadius  = tAllRadii(iIter);

            moris::ge::Sphere tSphere( tRadius, tXcenter, tYcenter, tZcenter );

            moris::ge::GEN_Phase_Table      tPhaseTable( 1,  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tSphere, tPhaseTable );

            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );
            for( uint i=0; i<tNumNodes; i++ )
            {
                tGENGeometryEngine.initialize_geometry_object_phase_values( tAllCoords );
                tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );
            }
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tFieldNames(iIter), EntityRank::NODE, tLSVals);
        }
        //------------------------------------------------------------------------------
//        std::string tOutputFile = "./sweep_radius.exo";
//        tInterpMesh->create_output_mesh(tOutputFile);

        delete tInterpMesh;
        delete tIntegMesh;
    }
}
//------------------------------------------------------------------------------
TEST_CASE("param_test_02","[GE],[param_sweep_02]")
{
    if(par_size()<=1)
    {
        size_t tModelDimension = 3;
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
        tParameters.set( "domain_dimensions",                "10, 10, 10" );
        tParameters.set( "domain_offset",                    "-5, -5, -5" );

//        tParameters.set( "domain_sidesets",      "1, 2, 3, 4, 5, 6" );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_output_meshes", "0" );
        tParameters.set( "lagrange_input_meshes", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );

        tParameters.insert( "initial_refinement", 0 );

        //  HMR Initialization
        moris::hmr::HMR tHMR( tParameters );

        std::shared_ptr< hmr::Database > tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement( 0 );

        tDatabase->update_bspline_meshes();
        tDatabase->update_lagrange_meshes();

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );
        tHMR.finalize();

        std::shared_ptr< hmr::Interpolation_Mesh_HMR >      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );

        mtk::Mesh_Manager tMesh1;
        //------------------------------------------------------------------------------
        real tXcenter = 0.0;
        real tYcenter = 0.0;
        real tZcenter = 0.0;

        Matrix< DDRMat > tAllRadii{{ 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                     2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9 }};
        uint tNumIters = tAllRadii.n_cols();
        //------------------------------------------------------------------------------
        for(uint iIter=0; iIter<tNumIters; iIter++)    // loop over radii values
        {
            real tRadius  = tAllRadii(iIter);
            moris::ge::Sphere tSphere( tRadius, tXcenter, tYcenter, tZcenter );

            moris::ge::GEN_Phase_Table      tPhaseTable( 1,  Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tSphere, tPhaseTable, tModelDimension );
            xtk::Model                      tXTKModel( tModelDimension, tInterpMesh.get(), tGENGeometryEngine );
            tXTKModel.mVerbose = false;

            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE_1,0);

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();

            // place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);

            xtk::Output_Options tOutputOptions;
            tOutputOptions.mAddNodeSets = false;
            tOutputOptions.mAddSideSets = false;
            tOutputOptions.mAddClusters = false;

            std::string tCircleFieldName = "circleLS";
            tOutputOptions.mRealNodeExternalFieldNames = { tCircleFieldName };

            moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh(tOutputOptions);
            //------------------------------------------------------------------------------
            // get level-set field and add to integration mesh output
            uint tNumNodes = tIntegMesh1->get_num_entities(EntityRank::NODE);

            moris::Cell< moris::mtk::Vertex const * > tAllVerices = tIntegMesh1->get_all_vertices();

            Matrix< DDRMat > tAllCoords( tNumNodes, 3 );
            for(uint i=0; i<tNumNodes; i++ )
            {
                tAllCoords.set_row( i,tAllVerices(i)->get_coords() );
            }

            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );
            tGENGeometryEngine.initialize_geometry_object_phase_values( tAllCoords );

            Matrix< DDRMat > tLSVals( tNumNodes, 1, -1.0 );
            for( uint i=0; i<tNumNodes; i++ )
            {
                tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );
            }
            tIntegMesh1->add_mesh_field_real_scalar_data_loc_inds( tCircleFieldName, EntityRank::NODE, tLSVals );

            //------------------------------------------------------------------------------
            // output solution to mesh
//            std::string s = std::to_string(iIter);
//            std::string tMeshOutputFile = "circle" + s +".e";
//
//            tIntegMesh1->create_output_mesh(tMeshOutputFile);
            //------------------------------------------------------------------------------
            delete tIntegMesh1;
            for(uint k=0; k<tAllVerices.size(); k++)
            {
                delete tAllVerices(k);
            }
        }
        //------------------------------------------------------------------------------
    }
}
}   // end ge namespace
}   // end moris namespace




