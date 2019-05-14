///*
// * cl_GE_test_analytic.cpp
// *
// *  Created on: Apr 10, 2019
// *      Author: sonne
// */
//
//#include "cl_GE_Core.hpp"
//#include "catch.hpp"
//
////------------------------------------------------------------------------------
//// GE includes
//#include "cl_GE_Factory.hpp"
//
//// linalg includes
//#include "cl_Matrix.hpp"
//#include "fn_all_true.hpp"
//#include "fn_equal_to.hpp"
//#include "linalg_typedefs.hpp"
//#include "op_equal_equal.hpp"
//
//// MTK includes
//#include "cl_MTK_Cell.hpp"
//#include "cl_MTK_Enums.hpp"
//#include "cl_Mesh_Factory.hpp"
//#include "cl_MTK_Mesh.hpp"
//#include "cl_MTK_Mesh_Data_Input.hpp"
//#include "cl_MTK_Mesh_Tools.hpp"
//#include "cl_MTK_Scalar_Field_Info.hpp"
//#include "cl_MTK_Vertex.hpp"
//
//// FEM includes
//#include "cl_FEM_Field_Interpolator.hpp"
//#include "cl_FEM_Integrator.hpp"
//
////------------------------------------------------------------------------------
//
//using namespace moris;
//using namespace ge;
//
//TEST_CASE("analytic_functionalities_test_2D","[GE],[analytic_functionalities_2D]")
//        {
//            /*
//             * 1) create a single-element mesh
//             * 2) create the Geometry Engine
//             *      2.1) create the field representation(s) (geometry objects)
//             *      2.2) set mesh and T-matrix for object(s)
//             *      ---------------- occurs when GE_Core is created ----------------
//             *      2.2.1) initialize reps, compute ADVs (L2 projection)
//             *      2.2.2) create output object (intersection object) for each field rep
//             *      --------------------------------------------------------
//             * 3) determine LS values and sensitivity at nodes from output object
//             * 4) determine intersection locations along edges from output object
//             *
//             *         [3]
//             *   (0,1)      (1,1)
//             *      x--------x
//             *      |        |
//             *      O        |
//             * [4]  |        |  [2]
//             *      |        |
//             *      x-----O--x
//             *   (0,0)      (1,0)
//             *         [1]
//             */
//
//            /*
//             * --------------------------------------------------------
//             * (1) create a single-element mesh
//             * --------------------------------------------------------
//             */
//            uint aNumElemTypes = 1;     // quad
//            uint aNumDim = 2;           // specify number of spatial dimensions
//
//            Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh
//
//            Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads
//
//            Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
//                                        { 1.0, 0.0 },
//                                        { 1.0, 1.0 },
//                                        { 0.0, 1.0 }};             // Node coordinate matrix
//
//            Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map
//            //------------------------------------------------------------------------------
//            // create MORIS mesh using MTK database
//            mtk::MtkMeshData tMeshData( aNumElemTypes );
//            tMeshData.CreateAllEdgesAndFaces = true;
//            tMeshData.SpatialDim = & aNumDim;
//            tMeshData.ElemConn(0) = & aElementConnQuad;
//            tMeshData.NodeCoords = & aCoords;
//            tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
//            tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
//            //------------------------------------------------------------------------------
//            // declare scalar node field for the circle LS
//            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//            std::string tFieldName = "circle";
//            tNodeField1.set_field_name(tFieldName);
//            tNodeField1.set_field_entity_rank(EntityRank::NODE);
//
//            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField2;
//            std::string tFieldName1 = "projectionVals";
//            tNodeField2.set_field_name(tFieldName1);
//            tNodeField2.set_field_entity_rank(EntityRank::NODE);
//
//            // initialize field information container
//            moris::mtk::MtkFieldsInfo tFieldsInfo;
//            // Place the node field into the field info container
//            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//            add_field_for_mesh_input(&tNodeField2,tFieldsInfo);
//
//            // declare some supplementary fields
//            tMeshData.FieldsInfo = &tFieldsInfo;
//
//            // create mesh pair
//            mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
//            mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);
//
//            // place the pair in mesh manager
//            mtk::Mesh_Manager* tMeshManager;
//            uint tMeshIndex = tMeshManager->register_mesh_pair(tInterpMesh1,tIntegMesh1);
//
//            Matrix< DDRMat > tTMatrix( 4,4, 0.0 ); //T-matrix to be used for L2 projection
//            tTMatrix(0,0) = 0.25;
//            tTMatrix(1,1) = 0.25;
//            tTMatrix(2,2) = 0.25;
//            tTMatrix(3,3) = 0.25;
//
//            // input parameters for the circle LS
//            moris::Cell< real > tCircleInputs(3);
//            tCircleInputs(0) = 0.0;   // x center
//            tCircleInputs(1) = 0.0;   // y center
//            tCircleInputs(2) = 0.6;   // radius
//            /*
//             * --------------------------------------------------------
//             * (2) create geometry engine
//             * --------------------------------------------------------
//             */
//            Ge_Factory tFactory; // geometry pointer and type, set LS function and sensitivity
//            std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(type::ANALYTIC);
//            tGeom1->set_analytical_function(type::CIRCLE);
//            tGeom1->set_analytical_function_dphi_dx(type::CIRCLE);
//
//            tGeom1->set_mesh_and_t_matrix(tMeshManager, tTMatrix);
//
//            GE_Core tGeometryEngine;
//            tGeometryEngine.set_geometry(tGeom1);
//
//            // determine LS values at the nodes
//            Matrix< DDRMat > tLSVals(4,1,0.0);
//            Cell< Matrix<DDRMat > > tSensitivities(4);
//
//            uint tNumNode = tMeshManager->get_integration_mesh(tMeshIndex)->get_num_entities(EntityRank::NODE);
//            for(uint i=0; i<tNumNode; i++)
//            {
//                tLSVals(i,0)      = tGeometryEngine.get_geometry_pointer(0)->get_field_val_at_coordinate( tMeshManager->get_integration_mesh(tMeshIndex)->get_node_coordinate(i), tCircleInputs );
//                tSensitivities(i) = tGeometryEngine.get_geometry_pointer(0)->get_sensitivity_dphi_dp_at_coordinate( tMeshManager->get_integration_mesh(tMeshIndex)->get_node_coordinate(i), tCircleInputs );
//            }
//            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);
//            /*
//             * --------------------------------------------------------
//             * check determined values
//             * --------------------------------------------------------
//             */
//            CHECK( equal_to( tLSVals( 0,0 ), -0.36 ) );
//            CHECK( equal_to( tLSVals( 1,0 ),  0.64 ) );
//            CHECK( equal_to( tLSVals( 2,0 ),  1.64 ) );
//            CHECK( equal_to( tLSVals( 3,0 ),  0.64 ) );
//
//            Matrix< DDRMat > tCheckMat(3,2);
//            tCheckMat(0,0) = 1.0;
//            tCheckMat(0,1) = 1.0;
//            tCheckMat(1,0) = 1.0;
//            tCheckMat(1,1) = 0.0;
//            tCheckMat(2,0) = 0.0;
//            tCheckMat(2,1) = 1.0;
//            bool tMatrixMatch = all_true( tSensitivities(0) == tCheckMat );
//            CHECK( tMatrixMatch );
//
//            REQUIRE( tSensitivities(1)(0,1) == Approx(-0.75));
//            REQUIRE( tSensitivities(2)(0,0) == Approx(-0.75));
//            REQUIRE( tSensitivities(2)(0,1) == Approx(-0.75));
//            REQUIRE( tSensitivities(3)(0,0) == Approx(-0.75));
//            /*
//             * --------------------------------------------------------
//             * determine intersection points
//             * --------------------------------------------------------
//             */
//            uint tWhichGeom = 0;
//
//            Matrix< DDRMat > tADVs = tGeometryEngine.compute_nodal_advs( tWhichGeom,
//                                                                         tCircleInputs,
//                                                                         mtk::Geometry_Type::QUAD,
//                                                                         fem::Integration_Type::GAUSS,
//                                                                         fem::Integration_Order::QUAD_2x2,
//                                                                         fem::Interpolation_Type::LAGRANGE,
//                                                                         mtk::Interpolation_Order::LINEAR );
//            Matrix< DDRMat > tProjectedVals = tTMatrix*tADVs;
//            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName1, EntityRank::NODE, tProjectedVals);
//            CHECK( equal_to( tProjectedVals( 0,0 ), -0.36 ) );
//            CHECK( equal_to( tProjectedVals( 1,0 ),  0.64 ) );
//            CHECK( equal_to( tProjectedVals( 2,0 ),  1.64 ) );
//            CHECK( equal_to( tProjectedVals( 3,0 ),  0.64 ) );
//
//            //------------------------------------------------------------------------------
////            std::string tOutputFile = "./analytic_functionalities_2D.exo";
////            tMesh2D_Quad4->create_output_mesh(tOutputFile);
//            //cleanup
////            delete tMesh2D_Quad4;
//
//        }
////------------------------------------------------------------------------------
//
//
//TEST_CASE("analytic_functionalities_test_3D","[GE],[analytic_functionalities_3D]")
//{
//    /*
//     * create 3D mtk mesh, add a sphere LS function, determine intersection locations
//     */
//    //------------------------------------------------------------------------------
//    const std::string tFileName2 = "generated:2x2x2"; // Define background mesh size
//
//    // Declare scalar node field 1
//    moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField1;
//    std::string tSphere1FieldName = "sphere1";
//    tNodeSphereField1.set_field_name(tSphere1FieldName);
//    tNodeSphereField1.set_field_entity_rank(EntityRank::NODE);
//
//    // Initialize field information container
//    moris::mtk::MtkFieldsInfo tFieldsInfo;
//
//    // Place the node field into the field info container
//    add_field_for_mesh_input(&tNodeSphereField1,tFieldsInfo);
//
//    // Declare some supplementary fields
//    mtk::MtkMeshData tMeshData;
//    tMeshData.FieldsInfo = &tFieldsInfo;
//
//    // Create MORIS mesh using MTK database
//    mtk::Mesh* tMesh3DHexs = mtk::create_mesh( MeshType::STK, tFileName2, & tMeshData );
//    //------------------------------------------------------------------------------
//
//    // Define parameters for sphere
//    moris::Cell< real > tInput1(4);
//    tInput1(0) = 1.00; // x location of center
//    tInput1(1) = 1.00; // y location of center
//    tInput1(2) = 1.00; // z location of center
//    tInput1(3) = 0.50; // radius of sphere
//
//    //------------------------------------------------------------------------------
//    Ge_Factory tFactory;
//    std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(type::ANALYTIC);
//    tGeom1->set_analytical_function( type::SPHERE );
//
//    // Compute nodal sphere values for mesh
//    uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
//    Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);
//
//    // Collect nodal sphere values
//    for(uint i=0; i<tNumNodes; i++)
//    {
//        Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);
//
//        tNodeSphere1Vals(i) = tGeom1->get_field_val_at_coordinate( tNodeCoord,tInput1 );
//    }
//    // add nodal sphere values to mesh
//    tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);
//
//    //------------------------------------------------------------------------------
//
//    std::shared_ptr<Geometry_Engine_Interface> tInterface = std::make_shared< GE_Core >();
//
//    tInterface->set_geometry( tGeom1 );
//
//    /* *********************************
//     * this test is not yet complete...
//     * need to determine the
//     * intersection location(s)
//     * *********************************
//     */
//
//
//    //------------------------------------------------------------------------------
////    std::string tOutputFile = "./analyticTest.exo";
////    tMesh3DHexs->create_output_mesh(tOutputFile);
//    //------------------------------------------------------------------------------
//    // cleanup
//    delete tMesh3DHexs;
//
//}
//
//
////------------------------------------------------------------------------------
//
