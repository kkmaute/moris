/*
 * cl_General_Functionalities.cpp
 *
 *  Created on: Oct 24, 2019
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

#include "cl_GEN_Circle.hpp"

// MTK includes ---------------------------------
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"

// XTK include ----------------------------------
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
real tSimpleFunction( Matrix< DDRMat > const & aCoeff,
                      Matrix< DDRMat > const & aParam )
{
    return aCoeff(0)*aParam(0)*aParam(0) - aCoeff(1)*aParam(1);
}


//-----------------------------------------------------------------------------
TEST_CASE("general_test_01","[GE],[field_evaluation]")
{
    if(par_size()<=1)
    {
    /*
     * create 4-element stk mesh
     */
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
    // create mesh pair
    mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
    mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
    //-----------------------------------------------------------------------------
    /*
     * create geometry engine, register mesh and field
     */
    Matrix< DDRMat >tCoeff{ {5, 2} };
    moris::ge::GEN_Field_User_Defined* tUserDefField = new moris::ge::GEN_Field_User_Defined( tSimpleFunction, tCoeff );

    moris::ge::GEN_Geometry_Engine tGeomEng;
    moris_index tFieldIndex  = tGeomEng.register_field( tUserDefField );
    moris_index tGEMeshIndex = tGeomEng.set_mesh( &tMeshManager );
    /*
     * evaluate and check field values at nodes
     */
    Matrix< DDRMat > tNodeVals;
    tGeomEng.calc_field_vals_at_nodes( tGEMeshIndex, tFieldIndex, tNodeVals );

    Matrix< DDRMat > aSolVec{ {0,5,20,18,3,1,-4,-2,16} };

    bool tMatrixMatch = all_true( tNodeVals == aSolVec );
    CHECK( tMatrixMatch );
    //-----------------------------------------------------------------------------
//    std::string tOutputFile = "./general_test.exo";
//    tInterpMesh1->create_output_mesh(tOutputFile);
    //-----------------------------------------------------------------------------
    // clean up
    delete tInterpMesh;
    delete tIntegMesh;

    delete tUserDefField;
    }
}

//-----------------------------------------------------------------------------
TEST_CASE("general_test_02","[GE],[property_class]")
{
    if(par_size()<=1)
    {
        /*
         * Steps for FD check:
         *  (1) at initial configuration, compute and store:
         *          xgamma
         *          dxgamma/dphi
         *
         *  (2) perturb up:
         *          phi_1 = phi_1_0 + epsilon
         *          store xgamma_u
         *
         *  (3) perturb down:
         *          phi_1 = phi_1_0 - epsilon
         *          store xgamma_d
         *
         *  (4) dxgamma/dphi (FD) = (xgamma_u - xgamma_d)/2*epsilon
         *
         *  (5) compare dxgamma/dphi to dxgamma/dphi (FD)
         */
        uint tNumElemTypes = 1;     // quad
        uint tNumDim = 2;           // specify number of spatial dimensions

        Matrix< IdMat > tElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh

        Matrix< IdMat > tElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads

        Matrix< DDRMat > tCoords = {{ 0.0, 0.0 },
                                    { 1.0, 0.0 },
                                    { 1.0, 1.0 },
                                    { 0.0, 1.0 }};             // Node coordinate matrix

        Matrix< IdMat > tNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map
        //------------------------------------------------------------------------------
        // create MORIS mesh using MTK database
        mtk::MtkMeshData tMeshData( tNumElemTypes );
        tMeshData.CreateAllEdgesAndFaces = true;
        tMeshData.SpatialDim = & tNumDim;
        tMeshData.ElemConn(0) = & tElementConnQuad;
        //------------------------------------------------------------------------------
        tMeshData.NodeCoords = & tCoords;
        tMeshData.LocaltoGlobalElemMap(0) = & tElemLocalToGlobalQuad;
        tMeshData.LocaltoGlobalNodeMap = & tNodeLocalToGlobal;
        //------------------------------------------------------------------------------
        // declare scalar node field for the circle LS
//        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
//        std::string tFieldName = "circle";
//        tNodeField1.set_field_name(tFieldName);
//        tNodeField1.set_field_entity_rank(EntityRank::NODE);
//
//        // initialize field information container
//        moris::mtk::MtkFieldsInfo tFieldsInfo;
//        // Place the node field into the field info container
//        add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
//
//        // declare some supplementary fields
//        tMeshData.FieldsInfo = &tFieldsInfo;

        // create mesh pair
        mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
        mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
        //------------------------------------------------------------------------------
        real tRadius  = 0.6;
        real tXcenter = 0;
        real tYcenter = 0;
        moris::ge::Circle tCircle( tRadius, tXcenter, tYcenter );

        moris::ge::GEN_Phase_Table      tPhaseTable( 1,  Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tCircle, tPhaseTable, tNumDim );

        //=================== manual calls to GE (w/out XTK model) =============================
        uint tNumNodes = tInterpMesh->get_num_nodes();
        tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );

        Matrix< DDRMat > tNodeCoord;
        for( uint i=0; i<tNumNodes; i++ )
        {
            tNodeCoord = tInterpMesh->get_node_coordinate( i );
            tGENGeometryEngine.initialize_geometry_object_phase_values( tNodeCoord );
        }

        moris::Matrix< moris::IndexMat > tNodetoEdgeConnectivity(4, 2);
        // Edge 0
        (tNodetoEdgeConnectivity)(0, 0) = 0;
        (tNodetoEdgeConnectivity)(0, 1) = 1;
        // Edge 1
        (tNodetoEdgeConnectivity)(1, 0) = 1;
        (tNodetoEdgeConnectivity)(1, 1) = 2;
        // Edge 2
        (tNodetoEdgeConnectivity)(2, 0) = 2;
        (tNodetoEdgeConnectivity)(2, 1) = 3;
        // Edge 3
        (tNodetoEdgeConnectivity)(3, 0) = 3;
        (tNodetoEdgeConnectivity)(3, 1) = 0;

        Cell< GEN_Geometry_Object > tGeometryObjects;
        tGENGeometryEngine.is_intersected( tCoords, tNodetoEdgeConnectivity, (size_t) 1, tGeometryObjects );

        size_t tNumIntersections = tGeometryObjects.size();
//        std::cout<<tNumIntersections<<std::endl;
        //======================================
//        CHECK( tNumIntersections == 2 );  // there should be two intersected edges
        //======================================

        //================================== end ===============================================

//        xtk::Model                      tXTKModel( tNumDim, tInterpMesh, tGENGeometryEngine );
//        tXTKModel.mVerbose = false;
//        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_HIERARCHY_TET4};
//        tXTKModel.decompose(tDecompositionMethods);
        //------------------------------------------------------------------------------

    }
}

}   // ge namespace
}   // moris namespace
