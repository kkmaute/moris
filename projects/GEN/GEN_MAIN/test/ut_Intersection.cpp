#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_XTK_Edge_Topology.hpp"

//#include "cl_HMR.hpp"
//#include "cl_HMR_Mesh.hpp"
//#include "cl_HMR_Mesh_Interpolation.hpp"
//#include "cl_HMR_Mesh_Integration.hpp"
//#include "cl_PRM_HMR_Parameters.hpp"
//#include "cl_MTK_Writer_Exodus.hpp"
//#include "cl_XTK_Model.hpp"
//#include "cl_XTK_Enriched_Integration_Mesh.hpp"
//#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

namespace moris
{
    namespace ge
    {
        TEST_CASE("Test for Intersection Coordinates as PDVs", "[GE], [INTERSECTION_2D]")
        {
            // Mesh parameters
            uint tNumElemTypes = 1;     // quad
            uint tNumDim = 2;           // specify number of spatial dimensions
            Matrix< IdMat > tElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh
            Matrix< IdMat > tElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads
            Matrix< DDRMat > tCoords = {{ 0.0, 0.0 }, // Node coordinates
                                        { 1.0, 0.0 },
                                        { 1.0, 1.0 },
                                        { 0.0, 1.0 }};
            Matrix< IdMat > tNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map

            // Create mesh using MTK database
            mtk::MtkMeshData tMeshData( tNumElemTypes );
            tMeshData.CreateAllEdgesAndFaces = true;
            tMeshData.SpatialDim = & tNumDim;
            tMeshData.ElemConn(0) = & tElementConnQuad;
            tMeshData.NodeCoords = & tCoords;
            tMeshData.LocaltoGlobalElemMap(0) = & tElemLocalToGlobalQuad;
            tMeshData.LocaltoGlobalNodeMap = & tNodeLocalToGlobal;

            // Declare scalar node field for the circle LS
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName = "circle";
            tNodeField1.set_field_name(tFieldName);
            tNodeField1.set_field_entity_rank(EntityRank::NODE);

            // Initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;

            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1, tFieldsInfo);

            // Declare some supplementary fields
            tMeshData.FieldsInfo = &tFieldsInfo;

            // Create mesh pair
            mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
            mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

            // Place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );

            // Create circle geometry
            real tRadius = 0.6;
            Matrix<DDRMat> tADVs = {{0.0, 0.0, tRadius}};
            Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
            tGeometry(0) = std::make_shared<moris::ge::Circle>(tADVs,
                                                               Matrix<DDUMat>({{0, 1, 2}}),
                                                               Matrix<DDUMat>({{0, 1, 2}}),
                                                               Matrix<DDRMat>(0, 0));

            // Create geometry engine
            moris::ge::Phase_Table      tPhaseTable( 1, moris::ge::Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::Geometry_Engine  tGeometryEngine( tGeometry, tPhaseTable, tNumDim );
            tGeometryEngine.register_mesh(&tMeshManager);

            //=================== manual calls to GE (w/out XTK model) =============================
            uint tNumNodes = tInterpMesh->get_num_nodes();
            tGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );

            Matrix< DDRMat > tLSVals( tNumNodes,1 );
            for( uint i=0; i<tNumNodes; i++ )
            {
                tGeometryEngine.initialize_geometry_object_phase_values( tCoords );
                tLSVals(i) = tGeometryEngine.get_entity_phase_val( i,0 );
            }
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);

            // Edge connectivity
            moris::Matrix< moris::IndexMat > tNodetoEdgeConnectivity(4, 2);
            (tNodetoEdgeConnectivity)(0, 0) = 0;
            (tNodetoEdgeConnectivity)(0, 1) = 1;
            (tNodetoEdgeConnectivity)(1, 0) = 1;
            (tNodetoEdgeConnectivity)(1, 1) = 2;
            (tNodetoEdgeConnectivity)(2, 0) = 2;
            (tNodetoEdgeConnectivity)(2, 1) = 3;
            (tNodetoEdgeConnectivity)(3, 0) = 3;
            (tNodetoEdgeConnectivity)(3, 1) = 0;

            // Check that there are 2 intersected edges
            Cell< GEN_Geometry_Object > tGeometryObjects;
            tGeometryEngine.is_intersected( tCoords, tNodetoEdgeConnectivity, (size_t) 1, tGeometryObjects );
            size_t tNumIntersections = tGeometryObjects.size();
            REQUIRE( tNumIntersections == 2 );

            Matrix<DDRMat> tLclCoord( tNumIntersections, 1 );
            Matrix<DDRMat> tNewNodeCoords( tNumIntersections,tNumDim );
            for (uint i = 0; i < tNumIntersections; i++)
            {
                tLclCoord(i) = tGeometryObjects(i).get_interface_lcl_coord();
                tNewNodeCoords.set_row( i,tGeometryObjects(i).get_interface_glb_coord() );
            }

            // Check coordinates
            CHECK(tLclCoord(0) == Approx(0.2));
            CHECK(tLclCoord(1) == Approx(-0.2));
            CHECK(tNewNodeCoords(0, 0) == Approx(tRadius));
            CHECK(tNewNodeCoords(0, 1) == Approx(0.0));
            CHECK(tNewNodeCoords(1, 0) == Approx(0.0));
            CHECK(tNewNodeCoords(1, 1) == Approx(tRadius));

            // Edge indices
            Matrix< moris::IndexMat > tEdgeIndices00{{ 0,1 }};
            std::shared_ptr< xtk::Topology > tTop00 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices00 );
            Matrix< moris::IndexMat > tEdgeIndices01{{ 0,3 }};
            std::shared_ptr< xtk::Topology > tTop01 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices01 );

            // XTK topology cell
            moris::Cell< xtk::Topology* > tTopCell(2);
            tTopCell(0) = tTop00.get();
            tTopCell(1) = tTop01.get();

            // // New node indices
            moris::Cell< moris_index > tIndices(2);
            tIndices(0) = 4;
            tIndices(1) = 5;

            // New node local coordinates
            moris::Cell< Matrix<DDRMat> > tLocalCoords(2);
            tLocalCoords(0) = {{0.2}};
            tLocalCoords(1) = {{-0.2}};

            // New node global coordinates
            moris::Cell< Matrix<DDRMat> > tGlbCoords(2);
            tGlbCoords(0) = tNewNodeCoords.get_row(0);
            tGlbCoords(1) = tNewNodeCoords.get_row(1);

            // Create geometry objects
            tGeometryEngine.create_new_node_geometry_objects( tIndices, true, tTopCell, tLocalCoords, tGlbCoords );

            // Create PDVs on integration mesh
            tGeometryEngine.create_ig_pdv_hosts();

            // Compute interface sensitivities
            Matrix< DDRMat > tAllCoords = {{ 0.0, 0.0 },
                                           { 1.0, 0.0 },
                                           { 1.0, 1.0 },
                                           { 0.0, 1.0 },
                                           { tRadius, 0.0 },
                                           { 0.0, tRadius }};
            Matrix< IndexMat > tInd{{ tIndices(0),tIndices(1) }};
            tGeometryEngine.compute_interface_sensitivity( tInd, tAllCoords, 0 );

            // Check sensitivities of geometry objects
            Matrix<DDRMat> tDxDadv = tGeometryEngine.get_node_dx_dp(4);
            REQUIRE(tDxDadv.n_rows() == 2);
            REQUIRE(tDxDadv.n_cols() == 3);
            CHECK(tDxDadv(0, 0) == Approx(tRadius));
            CHECK(tDxDadv(0, 1) == Approx(0.0));
            CHECK(tDxDadv(0, 2) == Approx(1.0));
            CHECK(tDxDadv(1, 0) == Approx(0.0));
            CHECK(tDxDadv(1, 1) == Approx(0.0));
            CHECK(tDxDadv(1, 2) == Approx(0.0));

            tDxDadv = tGeometryEngine.get_node_dx_dp(5);
            REQUIRE(tDxDadv.n_rows() == 2);
            REQUIRE(tDxDadv.n_cols() == 3);
            CHECK(tDxDadv(0, 0) == Approx(0.0));
            CHECK(tDxDadv(0, 1) == Approx(0.0));
            CHECK(tDxDadv(0, 2) == Approx(0.0));
            CHECK(tDxDadv(1, 0) == Approx(0.0));
            CHECK(tDxDadv(1, 1) == Approx(tRadius));
            CHECK(tDxDadv(1, 2) == Approx(1.0));

            // Clean up
            delete tInterpMesh;
            delete tIntegMesh;
        }

        //--------------------------------------------------------------------------------------------------------------

//        TEST_CASE("general_test_02","[GE],[sensitivity_check_02]")
//        {
//            if(par_size()<=1)
//            {
//                uint tLagrangeMeshIndex = 0;
//                std::string tFieldName = "Cylinder";
//
//                moris::ParameterList tParameters = prm::create_hmr_parameter_list();
//
//                tParameters.set( "number_of_elements_per_dimension", std::string("10, 10"));
//                tParameters.set( "domain_dimensions", std::string("2, 2") );
//                tParameters.set( "domain_offset", std::string("-1.0, -1.0") );
//                tParameters.set( "domain_sidesets", std::string("1,2,3,4") );
//                tParameters.set( "lagrange_output_meshes", std::string("0") );
//
//                tParameters.set( "lagrange_orders", std::string("1") );
//                tParameters.set( "lagrange_pattern", std::string("0") );
//                tParameters.set( "bspline_orders", std::string("1") );
//                tParameters.set( "bspline_pattern", std::string("0") );
//
//                tParameters.set( "lagrange_to_bspline", std::string("0") );
//
//                tParameters.set( "truncate_bsplines", 1 );
//                tParameters.set( "refinement_buffer", 3 );
//                tParameters.set( "staircase_buffer", 3 );
//                tParameters.set( "initial_refinement", 0 );
//
//                tParameters.set( "use_multigrid", 0 );
//                tParameters.set( "severity_level", 2 );
//
//                hmr::HMR tHMR( tParameters );
//
//                // initial refinement
//                tHMR.perform_initial_refinement( 0 );
//                tHMR.finalize();
//
//                moris::hmr::Interpolation_Mesh_HMR * tInterpolationMesh = tHMR.create_interpolation_mesh(tLagrangeMeshIndex);
//
//                Cell<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
//                tGeometry(0) = std::make_shared<moris::ge::Circle>(0.0, 0.0, 0.4501);
//
//                moris::ge::Phase_Table tPhaseTable (1, moris::ge::Phase_Table_Structure::EXP_BASE_2);
//                moris::ge::Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable, 2);
//
//                xtk::Model tXTKModel(2, tInterpolationMesh, &tGeometryEngine);
//                tXTKModel.mVerbose = false;
//
//                //Specify decomposition Method and Cut Mesh ---------------------------------------
//                Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//                tXTKModel.decompose(tDecompositionMethods);
//
//                tXTKModel.perform_basis_enrichment( EntityRank::NODE, 0 );
//                moris::print(tXTKModel.get_background_mesh().get_interface_nodes_loc_inds(0), "interface nodes");
//
//                xtk::Enriched_Interpolation_Mesh &tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
//                xtk::Enriched_Integration_Mesh &tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
//
//                // place the pair in mesh manager
//                mtk::Mesh_Manager tMeshManager;
//                tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);
//                tGeometryEngine.register_mesh(&tMeshManager);
//
//                // Calculate sensitivities on interface
//                Matrix<IndexMat> tInterfaceNodes = tXTKModel.get_background_mesh().get_interface_nodes_loc_inds(0);
//                Matrix<DDRMat> tCoordinates(tEnrIntegMesh.get_num_nodes(), 2);
//                for (uint tNodeIndex = 0; tNodeIndex < tEnrIntegMesh.get_num_nodes(); tNodeIndex++)
//                {
//                    tCoordinates.set_row(tNodeIndex, tEnrIntegMesh.get_node_coordinate(tNodeIndex));
//                }
//                tGeometryEngine.compute_interface_sensitivity( tInterfaceNodes,
//                                                               tCoordinates,
//                                                               0,
//                                                               false );
//
//                // Create PDVs on integration mesh
//                tGeometryEngine.create_ig_pdv_hosts();
//
//                // Make sure PDVs contain sensitivity information
//                Pdv_Host_Manager *tPdvHostManager = dynamic_cast<Pdv_Host_Manager *>(tGeometryEngine.get_design_variable_interface());
//                moris::print(tPdvHostManager->compute_diqi_dadv(), "sens");
//
//                delete tInterpolationMesh;
//
//            }
//        }
    }
}
