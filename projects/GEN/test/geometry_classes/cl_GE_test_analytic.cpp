/*
 * cl_GE_test_analytic.cpp
 *
 *  Created on: Apr 10, 2019
 *      Author: sonne
 */

#include "catch.hpp"


//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Main.hpp"
#include "cl_GE_Node.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

TEST_CASE("analytic_functionalities_test_2D","[GE],[analytic_functionalities_2D]")
        {
            /*
             * 1) create an element
             * 2) create analytic LS
             * 3) determine LS values and sensitivity at nodes
             * 4) determine intersection locations along edges
             *
             *         [3]
             *   (0,1)      (1,1)
             *      x--------x
             *      |        |
             *      O        |
             * [4]  |        |  [2]
             *      |        |
             *      x-----O--x
             *   (0,0)      (1,0)
             *         [1]
             */

            /*
             * --------------------------------------------------------
             * create GE node objects and use them to create GE element
             * --------------------------------------------------------
             */
//            mtk::Vertex* tVertex1 = new Node(0.0, 0.0);
//            mtk::Vertex* tVertex2 = new Node(1.0, 0.0);
//            mtk::Vertex* tVertex3 = new Node(1.0, 1.0);
//            mtk::Vertex* tVertex4 = new Node(0.0, 1.0);
//
//            moris::Cell< mtk::Vertex* > tNodes(4);
//            tNodes(0) = tVertex1;
//            tNodes(1) = tVertex2;
//            tNodes(2) = tVertex3;
//            tNodes(3) = tVertex4;
//
//            mtk::Cell* tElement = new Element(tNodes);
            uint aNumElemTypes = 1;     // quad
            uint aNumDim = 2;           // specify number of spatial dimensions

            Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh

            Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads

            Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
            { 1.0, 0.0 },
            { 1.0, 1.0 },
            { 0.0, 1.0 }};             // Node coordinate matrix

            Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map
            //------------------------------------------------------------------------------
            // create MORIS mesh using MTK database
            mtk::MtkMeshData tMeshData( aNumElemTypes );
            tMeshData.CreateAllEdgesAndFaces = true;
            tMeshData.SpatialDim = & aNumDim;
            tMeshData.ElemConn(0) = & aElementConnQuad;
            tMeshData.NodeCoords = & aCoords;
            tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
            tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
            //------------------------------------------------------------------------------
            mtk::Mesh* tMesh2D_Quad4 = create_mesh( MeshType::STK, tMeshData );
            /*
             * --------------------------------------------------------
             * create LS and determine values at nodes
             * --------------------------------------------------------
             */
            // geometry pointer and type, set LS function and sensitivity
            Ge_Factory tFactory;
            std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(type::ANALYTIC);
            tGeom1->set_analytical_function(type::CIRCLE);
            tGeom1->set_analytical_function_dphi_dx(type::CIRCLE);

            GE_Main tGeometryEngine;                    // create geometry engine and set geometry
            tGeometryEngine.set_geometry( tGeom1 );

            // input parameters for the circle
            moris::Cell< real > tCircleInputs(3);
            tCircleInputs(0) = 0.0;   // x center
            tCircleInputs(1) = 0.0;   // y center
            tCircleInputs(2) = 0.6;   // radius

            // determine LS values at the nodes
            Matrix< DDRMat > tLSVals(4,1,0.0);
            Cell< Matrix<DDRMat > > tSensitivities(4);

            uint tNumNode = tMesh2D_Quad4->get_num_entities(EntityRank::NODE);
            for(uint i=0; i<tNumNode; i++)
            {
                tLSVals(i,0)      = tGeometryEngine.get_geometry_pointer(0)->get_field_val_at_coordinate( tMesh2D_Quad4->get_node_coordinate(i), tCircleInputs );
                tSensitivities(i) = tGeometryEngine.get_geometry_pointer(0)->get_sensitivity_dphi_dp_at_coordinate( tMesh2D_Quad4->get_node_coordinate(i), tCircleInputs );
            }
            /*
             * --------------------------------------------------------
             * check determined values
             * --------------------------------------------------------
             */
            CHECK( equal_to( tLSVals( 0,0 ), -0.36 ) );
            CHECK( equal_to( tLSVals( 1,0 ),  0.64 ) );
            CHECK( equal_to( tLSVals( 2,0 ),  1.64 ) );
            CHECK( equal_to( tLSVals( 3,0 ),  0.64 ) );

            Matrix< DDRMat > tCheckMat(3,2);
            tCheckMat(0,0) = 1.0;
            tCheckMat(0,1) = 1.0;
            tCheckMat(1,0) = 1.0;
            tCheckMat(1,1) = 0.0;
            tCheckMat(2,0) = 0.0;
            tCheckMat(2,1) = 1.0;
            bool tMatrixMatch = all_true( tSensitivities(0) == tCheckMat );
            CHECK( tMatrixMatch );

            REQUIRE( tSensitivities(1)(0,1) == Approx(-0.75));
            REQUIRE( tSensitivities(2)(0,0) == Approx(-0.75));
            REQUIRE( tSensitivities(2)(0,1) == Approx(-0.75));
            REQUIRE( tSensitivities(3)(0,0) == Approx(-0.75));
            /*
             * --------------------------------------------------------
             * determine intersection points
             * --------------------------------------------------------
             */
            //T-matrix to be used for L2 projection
            Matrix< DDRMat > tTMat( 4,4, 0.0 );
            tTMat(0,0) = 0.25;
            tTMat(1,1) = 0.25;
            tTMat(2,2) = 0.25;
            tTMat(3,3) = 0.25;

            Matrix< DDRMat > tPhiHat;
std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;
            tPhiHat = tGeometryEngine.determine_phi_hat( tTMat, tMesh2D_Quad4, tCircleInputs, 0 );
print(tPhiHat,"phi hat");

            Matrix< DDRMat > tPhi = tTMat*tPhiHat;
            print(tPhi,"phi vals");

            //------------------------------------------------------------------------------
            //cleanup
//            delete tVertex1; delete tVertex2;
//            delete tVertex3; delete tVertex4;
//            delete tElement;
            delete tMesh2D_Quad4;

        }
//------------------------------------------------------------------------------


TEST_CASE("analytic_functionalities_test_3D","[GE],[analytic_functionalities_3D]")
{
    /*
     * create 3D mtk mesh, add a sphere LS function, determine intersection locations
     */
    //------------------------------------------------------------------------------
    const std::string tFileName2 = "generated:2x2x2"; // Define background mesh size

    // Declare scalar node field 1
    moris::mtk::Scalar_Field_Info<DDRMat> tNodeSphereField1;
    std::string tSphere1FieldName = "sphere1";
    tNodeSphereField1.set_field_name(tSphere1FieldName);
    tNodeSphereField1.set_field_entity_rank(EntityRank::NODE);

    // Initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;

    // Place the node field into the field info container
    add_field_for_mesh_input(&tNodeSphereField1,tFieldsInfo);

    // Declare some supplementary fields
    mtk::MtkMeshData tMeshData;
    tMeshData.FieldsInfo = &tFieldsInfo;

    // Create MORIS mesh using MTK database
    mtk::Mesh* tMesh3DHexs = mtk::create_mesh( MeshType::STK, tFileName2, & tMeshData );
    //------------------------------------------------------------------------------

    // Define parameters for sphere
    moris::Cell< real > tInput1(4);
    tInput1(0) = 1.00; // x location of center
    tInput1(1) = 1.00; // y location of center
    tInput1(2) = 1.00; // z location of center
    tInput1(3) = 0.50; // radius of sphere

    //------------------------------------------------------------------------------
    Ge_Factory tFactory;
    std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(type::ANALYTIC);
    tGeom1->set_analytical_function( type::SPHERE );

    // Compute nodal sphere values for mesh
    uint tNumNodes = tMesh3DHexs->get_num_entities(EntityRank::NODE);
    Matrix< DDRMat > tNodeSphere1Vals(1,tNumNodes);

    // Collect nodal sphere values
    for(uint i=0; i<tNumNodes; i++)
    {
        Matrix< DDRMat > tNodeCoord = tMesh3DHexs->get_node_coordinate(i);

        tNodeSphere1Vals(i) = tGeom1->get_field_val_at_coordinate( tNodeCoord,tInput1 );
    }
    // add nodal sphere values to mesh
    tMesh3DHexs->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tNodeSphere1Vals);

    //------------------------------------------------------------------------------

    std::shared_ptr<Geometry_Engine_Interface> tInterface = std::make_shared< GE_Main >();

    tInterface->set_geometry( tGeom1 );

    /* *********************************
     * this test is not yet complete...
     * need to determine the
     * intersection location
     * *********************************
     */


    //------------------------------------------------------------------------------

//    std::string tOutputFile = "./analyticTest.exo";
//    tMesh3DHexs->create_output_mesh(tOutputFile);
    //------------------------------------------------------------------------------
    // cleanup
    delete tMesh3DHexs;

}


//------------------------------------------------------------------------------

