///*
// * cl_General_Functionalities.cpp
// *
// *  Created on: Oct 24, 2019
// *      Author: sonne
// */
//#include "catch.hpp"
//
//#include "cl_Matrix.hpp"
//#include "fn_all_true.hpp"
//#include "op_equal_equal.hpp"
//// GE
//#include "../projects/GEN/src/new/additional/cl_GEN_Enums.hpp"
//#include "../projects/GEN/src/new/additional/cl_GEN_Phase_Table.hpp"
//
//#include "../projects/GEN/src/new/field/cl_GEN_Field_User_Defined.hpp"
//
//#include "../projects/GEN/src/new/geomeng/cl_GEN_Geometry_Engine.hpp"
//#include "../projects/GEN/src/new/geomeng/cl_GEN_Geometry_Engine.hpp"
//
//#include "../projects/GEN/src/new/geometry/cl_GEN_Circle.hpp"
//
//#include "cl_Mesh_Factory.hpp"
//#include "cl_MTK_Mesh_Manager.hpp"
//
////------------------------------------------------------------------------------
//
//namespace moris
//{
//namespace ge
//{
//real tSimpleFunction( Matrix< DDRMat > & aCoeff,
//                      Matrix< DDRMat > & aParam )
//{
//    return aCoeff(0)*aParam(0)*aParam(0) - aCoeff(1)*aParam(1);
//}
//
//
//
//TEST_CASE("general_tests","[GE],[field_evaluation]")
//{
//    if(par_size()<=1)
//    {
//    /*
//     * create 4-element stk mesh
//     */
////-----------------------------------------------------------------------------
//    uint aNumElemTypes = 1;     // quad
//    uint aNumDim = 2;           // specify number of spatial dimensions
//
//    Matrix< IdMat > aElementConnQuad = {{ 1, 2, 5, 8 },
//                                        { 2, 3, 4, 5 },
//                                        { 8, 5, 6, 7 },
//                                        { 5, 4, 9, 6 }};        // specify element connectivity of quad for mesh
//
//    Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1, 2, 3, 4 }};      // specify the local to global element map for quads
//
//    Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
//                                { 1.0, 0.0 },
//                                { 2.0, 0.0 },
//                                { 2.0, 1.0 },
//                                { 1.0, 1.0 },
//                                { 1.0, 2.0 },
//                                { 0.0, 2.0 },
//                                { 0.0, 1.0 },
//                                { 2.0, 2.0 }};      // Node coordinate matrix
//
//    Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};       // specify the local to global map
//    //------------------------------------------------------------------------------
//    // create MORIS mesh using MTK database
//    mtk::MtkMeshData tMeshData( aNumElemTypes );
//    tMeshData.CreateAllEdgesAndFaces = true;
//    tMeshData.SpatialDim = & aNumDim;
//    tMeshData.ElemConn(0) = & aElementConnQuad;
//    //------------------------------------------------------------------------------
//    tMeshData.NodeCoords = & aCoords;
//    tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
//    tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
//    //------------------------------------------------------------------------------
//    // declare scalar node field for the circle LS
//    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//    std::string tFieldName = "circle";
//    tNodeField1.set_field_name(tFieldName);
//    tNodeField1.set_field_entity_rank( EntityRank::NODE );
//
//    // initialize field information container
//    moris::mtk::MtkFieldsInfo tFieldsInfo;
//    // Place the node field into the field info container
//    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//
//    // declare some supplementary fields
//    tMeshData.FieldsInfo = &tFieldsInfo;
//
//    // create mesh pair
//    mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
//    mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );
//
//    // place the pair in mesh manager
//    mtk::Mesh_Manager tMeshManager;
//    uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
////-----------------------------------------------------------------------------
//    /*
//     * create geometry engine, register mesh and field
//     */
//    Matrix< DDRMat >tCoeff{ {5, 2} };
//    moris::ge::GEN_Field_User_Defined tUserDefField( tSimpleFunction, tCoeff );
//
//    moris::ge::GEN_Geometry_Engine tGeomEng;
//    moris_index tFieldIndex  = tGeomEng.register_field( &tUserDefField );
//    moris_index tGEMeshIndex = tGeomEng.set_mesh( &tMeshManager );
//
//    /*
//     * evaluate and check field values at the all nodes
//     */
//    Matrix< DDRMat > aNodeVals = tGeomEng.calc_field_vals_at_nodes( tGEMeshIndex, tFieldIndex );
//print( aNodeVals,"aNodeVals" );
//    Matrix< DDRMat > aSolVec{ {0,5,20,18,3,1,-4,-2,16} };
//print( aSolVec,"aSolVec" );
//    bool tMatrixMatch = all_true( aNodeVals == aSolVec );
//    CHECK( tMatrixMatch );
////-----------------------------------------------------------------------------
////    std::string tOutputFile = "./general_test.exo";
////    tInterpMesh1->create_output_mesh(tOutputFile);
//    // clean up
//    delete tInterpMesh;
//    delete tIntegMesh;
//    }
//}
//
//}   // ge namespace
//}   // moris namespace
//
//
