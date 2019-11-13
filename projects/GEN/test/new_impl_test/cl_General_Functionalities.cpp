/*
 * cl_General_Functionalities.cpp
 *
 *  Created on: Oct 24, 2019
 *      Author: sonne
 */

#include "catch.hpp"

#include "cl_Matrix.hpp"

// GE
#include "./ripped/additional/cl_GEN_Enums.hpp"
#include "./ripped/additional/cl_GEN_Phase_Table.hpp"

#include "./ripped/geomeng/cl_GEN_Geometry_Engine.hpp"

#include "./ripped/geometry/cl_GEN_Circle.hpp"

// MTK includes
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{

TEST_CASE("general_tests","[GE],[general_test_01]")
{
    if(par_size()<=1)
    {
    // create 4-element stk mesh
//-----------------------------------------------------------------------------
    uint aNumElemTypes = 1;     // quad
    uint aNumDim = 2;           // specify number of spatial dimensions

    Matrix< IdMat > aElementConnQuad = {{ 1, 2, 5, 8 },
                                        { 2, 3, 4, 5 },
                                        { 8, 5, 6, 7 },
                                        { 5, 4, 9, 6 }};        // specify element connectivity of quad for mesh

    Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1, 2, 3, 4 }};      // specify the local to global element map for quads

    Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
                                { 1.0, 0.0 },
                                { 2.0, 0.0 },
                                { 2.0, 1.0 },
                                { 1.0, 1.0 },
                                { 1.0, 2.0 },
                                { 0.0, 2.0 },
                                { 0.0, 1.0 },
                                { 2.0, 2.0 }};      // Node coordinate matrix

    Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8, 9 }};       // specify the local to global map
    //------------------------------------------------------------------------------
    // create MORIS mesh using MTK database
    mtk::MtkMeshData tMeshData( aNumElemTypes );
    tMeshData.CreateAllEdgesAndFaces = true;
    tMeshData.SpatialDim = & aNumDim;
    tMeshData.ElemConn(0) = & aElementConnQuad;
    //------------------------------------------------------------------------------
    tMeshData.NodeCoords = & aCoords;
    tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
    tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
    //------------------------------------------------------------------------------
    // declare scalar node field for the circle LS
    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
    std::string tFieldName = "circle";
    tNodeField1.set_field_name(tFieldName);
    tNodeField1.set_field_entity_rank( EntityRank::NODE );

    // initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    // Place the node field into the field info container
    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

    // declare some supplementary fields
    tMeshData.FieldsInfo = &tFieldsInfo;

    // create mesh pair
    mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
    mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh1 );

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);
//-----------------------------------------------------------------------------
    real tR = 0.6; real tX = 1.0; real tY = 1.0;

    ge::Circle tLSCircle( tR,tX,tY );

    GEN_Phase_Table tPhaseTable( 1, Phase_Table_Structure::EXP_BASE_2 );
    GEN_Geometry_Engine tGeometryEngine(tLSCircle,tPhaseTable);




//-----------------------------------------------------------------------------
//    std::string tOutputFile = "./general_test.exo";
//    tInterpMesh1->create_output_mesh(tOutputFile);
    // clean up
    delete tInterpMesh1;
    delete tIntegMesh1;
    }
}

}   // ge namespace
}   // moris namespace


