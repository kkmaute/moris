/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Model_2D.cpp
 *
 */

// #include "catch.hpp"
// #include "cl_Communication_Tools.hpp"

// #include "cl_XTK_Model.hpp"

// #include "cl_GEN_Circle.hpp"
// #include "cl_GEN_Geometry.hpp"

// using namespace moris;

// namespace xtk
// {

// TEST_CASE("2D Regular Subdivision Method","[RSM_2D_Lin]")
//                 {
//     int tProcSize = par_size();

//     if(tProcSize==1)
//     {
//         // Geometry Engine Setup -----------------------
//         // Using a Levelset Circle as the Geometry

//         real tRadius = 0.7;
//         real tXCenter = 1.0;
//         real tYCenter = 1.0;
//         Vector<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
//         tGeometry(0) = std::make_shared<moris::ge::Circle>(tXCenter, tYCenter, tRadius);

//         // Create Mesh ---------------------------------
//         // Generate data for test
//         uint aNumDim = 2;
//         Matrix< DDRMat >  aCoords(6,2);
//         aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
//         aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
//         aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
//         aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
//         aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
//         aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
//         Matrix< IdMat >     aElemConn( 2, 4 );

//         // 0D to 3D connectivity (node to element)
//         aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
//         aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

//         Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

//         // No need of an element map since elements in connectivity table are assumed to be contiguous

//         // Create MORIS mesh using MTK database
//         moris::mtk::MtkMeshData aMeshData;
//         aMeshData.CreateAllEdgesAndFaces = true;
//         aMeshData.SpatialDim = &aNumDim;
//         aMeshData.ElemConn(0)= &aElemConn;
//         aMeshData.NodeCoords = &aCoords;
//         aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//         // Setup XTK Model -----------------------------
//         size_t tModelDimension = 2;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         //Specify your decomposition methods and start cutting
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
//         tXTKModel.decompose(tDecompositionMethods);

//         // Access the decomposed XTK Mesh
//         Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

//         // Do some testing
//         size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
//         size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

//         CHECK(tNumNodesAfterDecompositionXTK == 10); /* two duplicates from the shared nodes*/
//         CHECK(tNumElementsAfterDecompositionXTK == 8);

//         moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//         moris::Matrix< moris::DDRMat > tExpectedNodeCoordinates = {{0.0, 0.0},
//                                                                    {1.0, 0.0},
//                                                                    {1.0, 1.0},
//                                                                    {0.0, 1.0},
//                                                                    {2.0, 0.0},
//                                                                    {2.0, 1.0},
//                                                                    {0.5, 0.5},
//                                                                    {1.5, 0.5}};
//         CHECK(equal_to(tNodeCoordinates,tExpectedNodeCoordinates));

//         moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
//         std::string tMeshOutputFile    = "./xtk_exo/xtk_test_output_regular_subdivision_2d.e";
//         tCutMeshData->create_output_mesh(tMeshOutputFile);
//         delete tMeshData;
//         delete tCutMeshData;

//     }
//                 }

// TEST_CASE("2D Conformal Subdivision","[CM_2D_LIN]")
// {
//     if(par_size()==1)
//     {
//         // Geometry Engine Setup -----------------------
//         // Using a Levelset Circle as the Geometry
//         real tRadius = 0.7;
//         real tXCenter = 1.0;
//         real tYCenter = 1.0;
//         Vector<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
//         tGeometry(0) = std::make_shared<moris::ge::Circle>(tXCenter, tYCenter, tRadius);
//         // Create Mesh ---------------------------------
//         // Generate data for test
//         uint aNumDim = 2;
//         Matrix< DDRMat >  aCoords(6,2);
//         aCoords(0,0) = 0.0, aCoords(0,1) = 0.0;
//         aCoords(1,0) = 1.0, aCoords(1,1) = 0.0;
//         aCoords(2,0) = 1.0, aCoords(2,1) = 1.0;
//         aCoords(3,0) = 0.0, aCoords(3,1) = 1.0;
//         aCoords(4,0) = 2.0, aCoords(4,1) = 0.0;
//         aCoords(5,0) = 2.0, aCoords(5,1) = 1.0;
//         Matrix< IdMat >     aElemConn( 2, 4 );

//         // 0D to 3D connectivity (node to element)
//         aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
//         aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 5; aElemConn( 1, 2 ) = 6; aElemConn( 1, 3 ) = 3;

//         Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

//         // No need of an element map since elements in connectivity table are assumed to be contiguous

//         // Create MORIS mesh using MTK database
//         moris::mtk::MtkMeshData aMeshData;
//         aMeshData.CreateAllEdgesAndFaces = true;
//         aMeshData.SpatialDim = &aNumDim;
//         aMeshData.ElemConn(0)= &aElemConn;
//         aMeshData.NodeCoords = &aCoords;
//         aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

// // Setup XTK Model ----------------------------------------------------------------
//         size_t tModelDimension = 2;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         //Specify decomposition Method and Cut Mesh ---------------------------------------
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//         tXTKModel.decompose(tDecompositionMethods);

//         // output to exodus file ----------------------------------------------------------
//         moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();

//         std::string tMeshOutputFile ="./xtk_exo/xtk_test_output_conformal_2d.e";
//         tCutMeshData->create_output_mesh(tMeshOutputFile);

//         delete tCutMeshData;
//         delete tMeshData;
//     }
// }

// TEST_CASE("2D Regular Subdivision Method Quadratic Lagrange Cells","[RSM_2D_Quad]")
//                 {
//     int tProcSize = par_size();

//     if(tProcSize==1)
//     {
//         // Geometry Engine Setup -----------------------
//         // Using a Levelset Circle as the Geometry

//         real tRadius = 0.7;
//         real tXCenter = 1.0;
//         real tYCenter = 1.0;
//         Vector<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
//         tGeometry(0) = std::make_shared<moris::ge::Circle>(tXCenter, tYCenter, tRadius);

//         // Create Mesh ---------------------------------
//         // Generate data for test
//         uint aNumDim = 2;
//         Matrix< DDRMat >  aCoords(15,2);
//         aCoords(0,0)  = 0.0, aCoords(0,1)  = 0.0;
//         aCoords(1,0)  = 1.0, aCoords(1,1)  = 0.0;
//         aCoords(2,0)  = 1.0, aCoords(2,1)  = 1.0;
//         aCoords(3,0)  = 0.0, aCoords(3,1)  = 1.0;
//         aCoords(4,0)  = 0.5, aCoords(4,1)  = 0.0;
//         aCoords(5,0)  = 1.0, aCoords(5,1)  = 0.5;
//         aCoords(6,0)  = 0.5, aCoords(6,1)  = 1.0;
//         aCoords(7,0)  = 0.0, aCoords(7,1)  = 0.5;
//         aCoords(8,0)  = 0.5, aCoords(8,1)  = 0.5;
//         aCoords(9,0)  = 2.0, aCoords(9,1)  = 0.0;
//         aCoords(10,0) = 2.0, aCoords(10,1) = 1.0;
//         aCoords(11,0) = 1.5, aCoords(11,1) = 0.0;
//         aCoords(12,0) = 2.0, aCoords(12,1) = 0.5;
//         aCoords(13,0) = 1.5, aCoords(13,1) = 1.0;
//         aCoords(14,0) = 1.5, aCoords(14,1) = 0.5;

//         // 0D to 3D connectivity (node to element)
//         Matrix< IdMat >     aElemConn( 2, 9 );
//         aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
//         aElemConn( 0, 4 ) = 5; aElemConn( 0, 5 ) = 6; aElemConn( 0, 6 ) = 7; aElemConn( 0, 7 ) = 8;
//         aElemConn( 0, 8 ) = 9;

//         aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 10 ; aElemConn( 1, 2 ) = 11; aElemConn( 1, 3 ) = 3;
//         aElemConn( 1, 4 ) = 12; aElemConn( 1, 5 ) = 13; aElemConn( 1, 6 ) = 14; aElemConn( 1, 7 ) = 6;
//         aElemConn( 1, 8 ) = 15;

//         Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

//         // Create MORIS mesh using MTK database
//         moris::mtk::MtkMeshData aMeshData;
//         aMeshData.CreateAllEdgesAndFaces = true;
//         aMeshData.SpatialDim = &aNumDim;
//         aMeshData.ElemConn(0)= &aElemConn;
//         aMeshData.NodeCoords = &aCoords;
//         aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//         std::string tMeshInterpOutputFile    = "./xtk_exo/xtk_test_output_regular_subdivision_interp_2d.e";
//         tMeshData->create_output_mesh(tMeshInterpOutputFile);

//         // Setup XTK Model -----------------------------
//         size_t tModelDimension = 2;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         //Specify your decomposition methods and start cutting
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
//         tXTKModel.decompose(tDecompositionMethods);

//         // Access the decomposed XTK Mesh
//         Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

//         // Do some testing
//         size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
//         size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

//         CHECK(tNumNodesAfterDecompositionXTK == 10); /* two duplicates from the shared nodes*/
//         CHECK(tNumElementsAfterDecompositionXTK == 8);

//         moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//         moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
//         std::string tMeshOutputFile    = "./xtk_exo/xtk_test_output_regular_subdivision_2d.e";
//         tCutMeshData->create_output_mesh(tMeshOutputFile);
//         delete tMeshData;
//         delete tCutMeshData;

//     }
//                 }

// TEST_CASE("2D Conformal Quadratic Lagrange Cells","[CM_2D_QUAD]")
//                 {
//     int tProcSize = par_size();

//     if(tProcSize==1)
//     {
//         // Geometry Engine Setup -----------------------
//         // Using a Levelset Circle as the Geometry

//         real tRadius = 0.7;
//         real tXCenter = 1.0;
//         real tYCenter = 1.0;
//         Vector<std::shared_ptr<moris::ge::Geometry>> tGeometry(1);
//         tGeometry(0) = std::make_shared<moris::ge::Circle>(tXCenter, tYCenter, tRadius);

//         // Create Mesh ---------------------------------
//         // Generate data for test
//         uint aNumDim = 2;
//         Matrix< DDRMat >  aCoords(15,2);
//         aCoords(0,0)  = 0.0, aCoords(0,1)  = 0.0;
//         aCoords(1,0)  = 1.0, aCoords(1,1)  = 0.0;
//         aCoords(2,0)  = 1.0, aCoords(2,1)  = 1.0;
//         aCoords(3,0)  = 0.0, aCoords(3,1)  = 1.0;
//         aCoords(4,0)  = 0.5, aCoords(4,1)  = 0.0;
//         aCoords(5,0)  = 1.0, aCoords(5,1)  = 0.5;
//         aCoords(6,0)  = 0.5, aCoords(6,1)  = 1.0;
//         aCoords(7,0)  = 0.0, aCoords(7,1)  = 0.5;
//         aCoords(8,0)  = 0.5, aCoords(8,1)  = 0.5;
//         aCoords(9,0)  = 2.0, aCoords(9,1)  = 0.0;
//         aCoords(10,0) = 2.0, aCoords(10,1) = 1.0;
//         aCoords(11,0) = 1.5, aCoords(11,1) = 0.0;
//         aCoords(12,0) = 2.0, aCoords(12,1) = 0.5;
//         aCoords(13,0) = 1.5, aCoords(13,1) = 1.0;
//         aCoords(14,0) = 1.5, aCoords(14,1) = 0.5;

//         // 0D to 3D connectivity (node to element)
//         Matrix< IdMat >     aElemConn( 2, 9 );
//         aElemConn( 0, 0 ) = 1; aElemConn( 0, 1 ) = 2; aElemConn( 0, 2 ) = 3; aElemConn( 0, 3 ) = 4;
//         aElemConn( 0, 4 ) = 5; aElemConn( 0, 5 ) = 6; aElemConn( 0, 6 ) = 7; aElemConn( 0, 7 ) = 8;
//         aElemConn( 0, 8 ) = 9;

//         aElemConn( 1, 0 ) = 2; aElemConn( 1, 1 ) = 10 ; aElemConn( 1, 2 ) = 11; aElemConn( 1, 3 ) = 3;
//         aElemConn( 1, 4 ) = 12; aElemConn( 1, 5 ) = 13; aElemConn( 1, 6 ) = 14; aElemConn( 1, 7 ) = 6;
//         aElemConn( 1, 8 ) = 15;

//         Matrix< IdMat >  aElemLocaltoGlobal = {{1},{2}};

//         // Create MORIS mesh using MTK database
//         moris::mtk::MtkMeshData aMeshData;
//         aMeshData.CreateAllEdgesAndFaces = true;
//         aMeshData.SpatialDim = &aNumDim;
//         aMeshData.ElemConn(0)= &aElemConn;
//         aMeshData.NodeCoords = &aCoords;
//         aMeshData.LocaltoGlobalElemMap(0) = &aElemLocaltoGlobal;

//         moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

//         moris::ge::Geometry_Engine_Parameters tGeometryEngineParameters;
//         tGeometryEngineParameters.mGeometries = tGeometry;
//         moris::ge::Geometry_Engine tGeometryEngine(tMeshData, tGeometryEngineParameters);

//         std::string tMeshInterpOutputFile    = "./xtk_exo/xtk_test_conformal_subdivision_quad_interp_2d.e";
//         tMeshData->create_output_mesh(tMeshInterpOutputFile);

//         // Setup XTK Model -----------------------------
//         size_t tModelDimension = 2;
//         Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//         tXTKModel.mVerbose  =  false;

//         //Specify your decomposition methods and start cutting
//         Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3};
//         tXTKModel.decompose(tDecompositionMethods);

//         // Access the decomposed XTK Mesh
//         Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();

//         // Do some testing
//         size_t tNumNodesAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::NODE);
//         size_t tNumElementsAfterDecompositionXTK = tCutMesh.get_num_entities(EntityRank::ELEMENT);

//         CHECK(tNumNodesAfterDecompositionXTK == 16);
//         CHECK(tNumElementsAfterDecompositionXTK == 16);

//         moris::Matrix< moris::DDRMat > tNodeCoordinates = tXTKModel.get_background_mesh().get_all_node_coordinates_loc_inds();

//         moris::mtk::Mesh* tCutMeshData = tXTKModel.get_output_mesh();
//         std::string tMeshOutputFile    = "./xtk_exo/xtk_test_output_conformal_quad_2d.e";
//         tCutMeshData->create_output_mesh(tMeshOutputFile);
//         delete tMeshData;
//         delete tCutMeshData;

//     }
//                 }
// }

