/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Ghost_Stabilization.cpp
 *
 */

// #include "catch.hpp"
//
// #include "cl_XTK_Model.hpp"
//
// #include "cl_XTK_Enriched_Integration_Mesh.hpp"
// #include "cl_MTK_Writer_Exodus.hpp"
//
//
// #include "cl_GEN_Geometry.hpp"
// #include "cl_GEN_Plane.hpp"
//
// #include "cl_MTK_Mesh_Factory.hpp"
//
// namespace moris::xtk
//{
//
// TEST_CASE("Face oriented ghost stabilization","[GHOST]")
//{
////    if(par_size() == 1)
////    {
//    moris::Matrix<DDRMat> tCenters = {{ 2.0,2.0,2.1 }};
//    moris::Matrix<DDRMat> tNormals = {{ 1.0,1.0,1.0 }};
//
//    Vector<std::shared_ptr<moris::gen::Level_Set_Geometry>> tGeometry(1);
//    tGeometry(0) = std::make_shared<moris::gen::Plane>(tCenters(0), tCenters(1), tCenters(2), tNormals(0), tNormals(1), tNormals(2));
//
//    moris::gen::Phase_Table tPhaseTable (tGeometry.size());
//    moris::gen::Geometry_Engine tGeometryEngine(tGeometry, tPhaseTable);
//
//    // Create Mesh ---------------------------------
//    std::string tMeshFileName = "generated:1x1x16|sideset:z";
//    moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshFileName );
//
//    // create model
//    size_t tModelDimension = 3;
//    Model tXTKModel(tModelDimension,tMeshData,&tGeometryEngine);
//    tXTKModel.mVerbose  =  true;
//
//    // decompose
//    Vector<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//
//    tXTKModel.perform_basis_enrichment( mtk::EntityRank::NODE,0);
//
//    tXTKModel.construct_face_oriented_ghost_penalization_cells();
//
//    tXTKModel.get_ghost_stabilization().visualize_ghost_on_mesh(0);
//
//    // Write mesh
//    moris::mtk::Writer_Exodus writer(&tXTKModel.get_enriched_integ_mesh());
//    writer.write_mesh("", "xtk_test_ghost.exo","","xtk_temp.exo");
//
//    // Write the fields
//    writer.set_time(0.0);
//    writer.close_file();
//
//
//    delete tMeshData;
////    }
//
//
//}
//}
