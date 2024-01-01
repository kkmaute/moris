/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Mesh.cpp
 *
 */

// #include <memory>   // Shared ptrs

// #include "catch.hpp"

// // XTKL: Logging and Assertion Includes
// #include "cl_Logger.hpp"

// // XTKL: Container includes
// #include "cl_Vector.hpp"

// // XTKL: Linear Algebra Includes

// #include "cl_Matrix.hpp"

// #include "catch.hpp"
// #include "cl_Matrix.hpp"
// #include "linalg_typedefs.hpp"
// #include "fn_tet_volume.hpp"
// #include "fn_equal_to.hpp"

// //XTK  Includes
// #include "cl_XTK_Cut_Mesh.hpp"

// #include "cl_XTK_Child_Mesh.hpp"
// #include "cl_XTK_Model.hpp"
// #include "cl_XTK_Output_Options.hpp"
// #include "fn_compute_xtk_model_volumes.hpp"
// #include "Child_Mesh_Verification_Utilities.hpp"

// #include "cl_GEN_Sphere.hpp"

// namespace xtk
// {

//     TEST_CASE("Simple Mesh Testing","[XTK][CUT_MESH]"){
//         // Functionality tested in this case
//         // - get_num_simple_meshes
//         // - set/get_node_inds/ids
//         // TODO Add add_index function test

//         //       Indices
//         //    9----10----11
//         //   /|    /|    /|    // Element 1 Indices - 0,1,4,3,6,7,10,9
//         //  6-----7-----8 |    // Element 2 Indices - 1,2,5,4,7,8,11,10
//         //  | |   | |   | |
//         //  | 3---|-4---|-5
//         //  |/    |/    |/
//         //  0-----1-----2
//         //

//         //         Ids
//         //    2----10----11
//         //   /|    /|    /|
//         //  7-----36----8 |    // Element 1 Ids - 14,5,18,4,7,36,10,2
//         //  | |   | |   | |    // Element 2 Ids - 5,3,8,18,36,8,11,10
//         //  | 4---|18---|-8
//         //  |/    |/    |/
//         //  14----5-----3

//         // Initialize an Cut Mesh with 2 simple meshes that are 3d
//         size_t tModelDim = 3;
//         size_t tNumSimpleMesh = 2;
//         Cut_Mesh tCutMesh(nullptr,tNumSimpleMesh, tModelDim);

//         // Make sure the correct number of simple meshes have been created
//         REQUIRE(tCutMesh.get_num_child_meshes() == 2);

//         // Add node Indices then node ids for each element
//         Child_Mesh & tCM1 = tCutMesh.get_child_mesh(0);// Index of element 1
//         Child_Mesh & tCM2 = tCutMesh.get_child_mesh(1);// Index of element 2
//         moris::Matrix< moris::IndexMat > tInds1({{0, 1, 4, 3, 6, 7, 10, 9}}); // Indices of element 1
//         moris::Matrix< moris::IdMat > tIds1({{14, 5, 18, 4, 7, 36, 10, 2}}); // Ids of element 1
//         moris::Matrix< moris::IndexMat > tInds2({{1, 2, 5, 4, 7, 8, 11, 10}}); // Indices of element 2
//         moris::Matrix< moris::IdMat > tIds2({{   5, 3, 8, 18, 36, 8, 11, 10}}); // Ids of element 2

//         // Set node indices
//         tCM1.add_node_indices(tInds1);
//         tCM2.add_node_indices(tInds2);

//         // Set node Ids
//         tCM1.add_node_ids(tInds1);
//         tCM2.add_node_ids(tInds2);

//         //Test node indices for element 1
//         moris::Matrix< moris::IndexMat > const & tNodeInds1 = tCM1.get_node_indices();
//         moris::Matrix< moris::IndexMat > const & tNodeInds2 = tCM2.get_node_indices();
//         moris::Matrix< moris::IdMat > const & tNodeIds1 = tCM1.get_node_ids();
//         moris::Matrix< moris::IdMat > const & tNodeIds2 = tCM2.get_node_ids();

//         REQUIRE(tNodeInds1(0, 0) == Approx(0));
//         REQUIRE(tNodeInds1(0, 1) == Approx(1));
//         REQUIRE(tNodeInds1(0, 2) == Approx(4));
//         REQUIRE(tNodeInds1(0, 3) == Approx(3));
//         REQUIRE(tNodeInds1(0, 4) == Approx(6));
//         REQUIRE(tNodeInds1(0, 5) == Approx(7));
//         REQUIRE(tNodeInds1(0, 6) == Approx(10));
//         REQUIRE(tNodeInds1(0, 7) == Approx(9));

//         //Test node indices for element 2
//         REQUIRE(tNodeInds2(0, 0) == 1);
//         REQUIRE(tNodeInds2(0, 1) == 2);
//         REQUIRE(tNodeInds2(0, 2) == 5);
//         REQUIRE(tNodeInds2(0, 3) == 4);
//         REQUIRE(tNodeInds2(0, 4) == 7);
//         REQUIRE(tNodeInds2(0, 5) == 8);
//         REQUIRE(tNodeInds2(0, 6) == 11);
//         REQUIRE(tNodeInds2(0, 7) == 10);

//         //Test node ids for element 1
//         REQUIRE(tNodeIds1(0, 0) == 0);
//         REQUIRE(tNodeIds1(0, 1) == 1);
//         REQUIRE(tNodeIds1(0, 2) == 4);
//         REQUIRE(tNodeIds1(0, 3) == 3);
//         REQUIRE(tNodeIds1(0, 4) == 6);
//         REQUIRE(tNodeIds1(0, 5) == 7);
//         REQUIRE(tNodeIds1(0, 6) == 10);
//         REQUIRE(tNodeIds1(0, 7) == 9);

//         //Test node ids for element 2
//         REQUIRE(tNodeIds2(0, 0) == 1);
//         REQUIRE(tNodeIds2(0, 1) == 2);
//         REQUIRE(tNodeIds2(0, 2) == 5);
//         REQUIRE(tNodeIds2(0, 3) == 4);
//         REQUIRE(tNodeIds2(0, 4) == 7);
//         REQUIRE(tNodeIds2(0, 5) == 8);
//         REQUIRE(tNodeIds2(0, 6) == 11);
//         REQUIRE(tNodeIds2(0, 7) == 10);

//     }

//     /*
//      * Tests the following:
//      * - Signed volume of regular subdivision
//      *
//      * Will test: signed surface area to make sure it's 0
//      */
//     TEST_CASE("Regular Subdivision Geometry Check","[VOLUME_CHECK_REG_SUB]")
//     {
//         // Geometry Engine Setup -----------------------
//         moris::real tRadius  = 0.5;
//         moris::real tXCenter = 1.0;
//         moris::real tYCenter = 1.0;
//         moris::real tZCenter = 2.0;

//         // Create Mesh ---------------------------------
//         std::string tMeshFileName = "generated:1x1x4";
//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName );
//         std::string tBackgroundMeshOutput = "./xtk_exo/volume_check_rs_bm.e";
//         tMeshData->create_output_mesh(tBackgroundMeshOutput);

//         Vector< std::shared_ptr<moris::ge::Geometry> > tGeometryVector(1);
//         tGeometryVector(0) = std::make_shared<moris::ge::Sphere>(tXCenter, tYCenter, tZCenter, tRadius);

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometryVector;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//         // Setup XTK Model -----------------------------
//         size_t tModelDimension = 3;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         // Specify your decomposition methods and start cutting
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
//         tXTKModel.decompose(tDecompositionMethods);
//         moris::real tGoldVolume = 4;

//         xtk::Output_Options tOutputOptions;
//         tOutputOptions.mAddParallelFields = true;

//         verify_child_mesh_ancestry(tXTKModel.get_background_mesh(),
//                 tXTKModel.get_cut_mesh());

//         //
//         moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//         // Compute volume
//         moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//         moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
//         moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//         moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);

//         moris::real tMyVolume = tParentPhase0Vol + tParentPhase1Vol + tChildPhase0Vol + tChildPhase1Vol;

//         // Collect all volumes
//         moris::real tGlbVolume = sum_all(tMyVolume);

//         CHECK(moris::equal_to(tGoldVolume,tGlbVolume));

//         /*
//          * Check surface area
//          */

//         delete tMeshData;
//     }

//     TEST_CASE("Node Hierarchy Volume Check","[VOLUME_CHECK_NH]")
//     {
//         // Geometry Engine Setup -----------------------
//         moris::real tRadius  = 3.1;
//         moris::real tXCenter = 1.0;
//         moris::real tYCenter = 1.0;
//         moris::real tZCenter = 6.0;

//         // Create Mesh ---------------------------------
//         std::string tMeshFileName = "generated:1x1x4";
//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, NULL );

//         Vector< std::shared_ptr<moris::ge::Geometry> > tGeometryVector(1);
//         tGeometryVector(0) = std::make_shared<moris::ge::Sphere>(tXCenter, tYCenter, tZCenter, tRadius);

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometryVector;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//         // Setup XTK Model -----------------------------
//         size_t tModelDimension = 3;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         // Specify your decomposition methods and start cutting
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//         tXTKModel.decompose(tDecompositionMethods);

//         verify_child_mesh_ancestry(tXTKModel.get_background_mesh(),
//                 tXTKModel.get_cut_mesh());

//         //
//         moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//         // Compute volume
//         moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//         moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
//         moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//         moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);

//         moris::real tMyVolume = tParentPhase0Vol + tParentPhase1Vol + tChildPhase0Vol + tChildPhase1Vol;

//         // Collect all volumes
//         moris::real tGoldVolume = 4;
//         moris::real tGlbVolume = sum_all(tMyVolume);

//         CHECK(moris::equal_to(tGoldVolume,tGlbVolume));

//         /*
//          * Check surface area
//          */
//         delete tMeshData;

//     }

//     TEST_CASE("Node Hierarchy Geometry Check","[VOLUME_CHECK_REG_SUB]")
//     {
//         if(par_size() == 1)
//         {
//             // Geometry Engine Setup -----------------------
//             // Using a Levelset Sphere as the Geometry
//             moris::real tRadius = 0.25;
//             moris::real tXCenter = 1.0;
//             moris::real tYCenter = 1.0;
//             moris::real tZCenter = 0;

//             /*
//              * Specify Mesh Inputs
//              */
//             // Create Mesh ---------------------------------
//             std::string tMeshFileName = "generated:1x1x1";
//             moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName, NULL );

//             Vector< std::shared_ptr<moris::ge::Geometry> > tGeometryVector(1);
//             tGeometryVector(0) = std::make_shared<moris::ge::Sphere>(tXCenter, tYCenter, tZCenter, tRadius);

//             moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//             tGeometryEngineParameters.mGeometries = tGeometryVector;
//             moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//             // Setup XTK Model -----------------------------
//             size_t tModelDimension = 3;
//             Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);

//             //Specify your decomposition methods and start cutting
//             Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
//             tXTKModel.decompose(tDecompositionMethods);

//             moris::Matrix<moris::DDRMat> tNodeCoords = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//             moris::real tParentPhase0Vol = compute_non_intersected_parent_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//             moris::real tParentPhase1Vol = compute_non_intersected_parent_element_volume_by_phase(1,tNodeCoords,tXTKModel);
//             moris::real tChildPhase0Vol  = compute_child_element_volume_by_phase(0,tNodeCoords,tXTKModel);
//             moris::real tChildPhase1Vol  = compute_child_element_volume_by_phase(1,tNodeCoords,tXTKModel);

//             moris::real tGoldVolume = 1;
//             CHECK(tGoldVolume == Approx(tParentPhase0Vol + tParentPhase1Vol + tChildPhase0Vol + tChildPhase1Vol));

//             delete tMeshData;
//         }

//     }

// }

