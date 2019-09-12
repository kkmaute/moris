/*
 * cl_GE_test_analytic.cpp
 *
 *  Created on: Apr 10, 2019
 *      Author: sonne
 */

#include "catch.hpp"

//------------------------------------------------------------------------------
// FEM includes
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Interpolation_Rule.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"
#include "cl_GE_Node.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Vertex.hpp"

//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

TEST_CASE("analytic_functionalities_test_2D","[GE],[analytic_functionalities_2D]")
        {
            /*
             * 1) create a single-element mesh
             * 2) create the Geometry Engine
             * 3) determine LS values and sensitivity at nodes
             *    3.1) add a vertex and check values
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
             * (1) create a single-element mesh
             * --------------------------------------------------------
             */
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
            //------------------------------------------------------------------------------
            tMeshData.NodeCoords = & aCoords;
            tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
            tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
            //------------------------------------------------------------------------------
            // declare scalar node field for the circle LS
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
            std::string tFieldName = "circle";
            tNodeField1.set_field_name(tFieldName);
            tNodeField1.set_field_entity_rank(EntityRank::NODE);

            // initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;
            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

            // declare some supplementary fields
            tMeshData.FieldsInfo = &tFieldsInfo;

            // create mesh pair
            mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
            mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);

            // place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);

            /*
             * --------------------------------------------------------
             * (2) create geometry engine
             * --------------------------------------------------------
             */
            // input parameters for the circle LS
            moris::Cell< real > tCircleInputs = {{0},{0},{0.6}};
            //------------------------------------------------------------------------------

            Ge_Factory tFactory;
            std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(GeomType::ANALYTIC);

            tGeom1->set_my_mesh(&tMeshManager);
            tGeom1->set_my_constants(tCircleInputs);

            tGeom1->set_analytical_function(AnalyticType::CIRCLE);
            tGeom1->set_analytical_function_dphi_dp(AnalyticType::CIRCLE);

            GE_Core tGeometryEngine;
            moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );
            /*
             * --------------------------------------------------------
             * (3) access LS values and sensitivity at nodes from the pdv info object
             * --------------------------------------------------------
             */
            PDV_Info* tPDVInfo = tGeometryEngine.get_pdv_info_pointer( tMyGeomIndex );

            Matrix< DDRMat > tLSVals                = tPDVInfo->get_field_vals(  );           // phi
            Cell< Matrix< DDRMat > > tSensitivities = tPDVInfo->get_sensitivity_vals(  );     // dphi/dp

            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);      // add the determined values as a field on the mesh (for output purposes)
            /*
             * --------------------------------------------------------
             * 3.1) add a node/value and re-check
             * --------------------------------------------------------
             */
            Node tNewNode(0.5,0.5);
            tNewNode.set_index( 29 );

            tPDVInfo->add_vertex_and_value( tNewNode );

            Matrix< DDRMat > tNodeVal = tGeometryEngine.get_field_vals( tMyGeomIndex,tNewNode.get_index() );
            REQUIRE( tNodeVal(0,0) == Approx( 0.1071 ) );

            //------------------------------------------------------------------------------
            /*
             * ------------------------------------------------------------
             * 4) determine the intersection points along edges [1] and [4]
             * ------------------------------------------------------------
             */

            // edge [1]:
            Matrix< DDRMat > tGlobalPos = {{0,0},{1,0}};    // edge [1] goes form x=0 to x=1
            Matrix< DDRMat > tTHat = {{0},{1}};
//            Matrix< DDRMat > tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(1)(0,0) }};
            Matrix< DDRMat > tUHat = {{ tLSVals(0) },{ tLSVals(1) }};

            Intersection_Object_Line tIntersectionObject;

            tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
            moris_index tXInd = tPDVInfo->compute_intersection( &tIntersectionObject );

            Matrix< F31RMat > tIntersectionAlongX = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd );

            // edge [4]:
            tGlobalPos = {{0,0},{0,1}};     // edge [2] goes from y=0 to y=1
            tTHat = {{0},{1}};
//            tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(3)(0,0) }};
            tUHat = {{ tLSVals(0) },{ tLSVals(3) }};

            tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
            moris_index tYInd = tPDVInfo->compute_intersection( &tIntersectionObject );

            Matrix< F31RMat > tIntersectionAlongY = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tYInd );
            //------------------------------------------------------------------------------
            /*
             * --------------------------------------------------------
             * check determined values
             * --------------------------------------------------------
             */
//            CHECK( equal_to( tLSVals(0)(0,0), -0.60 ) );
//            CHECK( equal_to( tLSVals(1)(0,0),  0.40 ) );
//            CHECK( equal_to( tLSVals(2)(0,0),  (std::sqrt(2)-0.6)) );
//            CHECK( equal_to( tLSVals(3)(0,0),  0.40 ) );
            CHECK( equal_to( tLSVals(0), -0.60 ) );
            CHECK( equal_to( tLSVals(1),  0.40 ) );
            CHECK( equal_to( tLSVals(2),  (std::sqrt(2)-0.6)) );
            CHECK( equal_to( tLSVals(3),  0.40 ) );

            Matrix< DDRMat > tCheckMat(3,2);
            tCheckMat(0,0) = 1.0;            tCheckMat(0,1) = 1.0;
            tCheckMat(1,0) = 1.0;            tCheckMat(1,1) = 0.0;
            tCheckMat(2,0) = 0.0;            tCheckMat(2,1) = 1.0;

            bool tMatrixMatch = all_true( tSensitivities(0) == tCheckMat );
            CHECK( tMatrixMatch );

            REQUIRE( tSensitivities(1)(0,1) == Approx(-0.75));
            REQUIRE( tSensitivities(2)(0,0) == Approx(-0.75));
            REQUIRE( tSensitivities(2)(0,1) == Approx(-0.75));
            REQUIRE( tSensitivities(3)(0,0) == Approx(-0.75));

            CHECK( equal_to( tIntersectionAlongX(0,0), 0.6 ) );
            CHECK( equal_to( tIntersectionAlongY(0,1), 0.6 ) );
            //------------------------------------------------------------------------------

            tPDVInfo->compute_intersection_sensitivity( &tIntersectionObject, tXInd );

            Matrix< DDRMat > tTemp = tPDVInfo->get_intersection_sensitivity( &tIntersectionObject );

//            Matrix< DDRMat > tTemp2 = tPDVInfo->get_dxgamma_dp( &tIntersectionObject );

            print( tTemp, "intersection sensitivities w.r.t the field, at x-intersection point" );
//            print( tTemp2, "intersection sensitivity w.r.t pdv, at x-intersection point " );

//            std::string tOutputFile = "./analytic_functionalities_2D.exo";
//            tInterpMesh1->create_output_mesh(tOutputFile);
            // clean up
            delete tInterpMesh1;
            delete tIntegMesh1;

        }
//------------------------------------------------------------------------------
