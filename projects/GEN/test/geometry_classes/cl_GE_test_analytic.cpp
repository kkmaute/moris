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
            tGeom1->set_analytical_function_dphi_dx(AnalyticType::CIRCLE);

            GE_Core tGeometryEngine;
            moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );
            /*
             * --------------------------------------------------------
             * (3) determine LS values and sensitivity at nodes
             * --------------------------------------------------------
             */
            Cell< Matrix< DDRMat > > tLSVals(4);
            Cell< Matrix< DDRMat > > tSensitivities(4);

            Matrix< DDRMat > tFieldData(4,1, 0.0);
            for(moris_index n=0; n<4; n++)
            {
                tLSVals(n)        = tGeometryEngine.get_field_vals(tMyGeomIndex,n);
                tSensitivities(n) = tGeometryEngine.get_sensitivity_vals(tMyGeomIndex,n);

                tFieldData(n) = tLSVals(n)(0,0);
            }
            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tFieldData);   // add the determined values as a field on the mesh
            /*
             * --------------------------------------------------------
             * 3.1) add a node/value and re-check
             * --------------------------------------------------------
             */
            Node tNewNode(0.5,0.5);
            tNewNode.set_index( 29 );

            tGeometryEngine.add_vertex_and_value( tNewNode, tMyGeomIndex );

            Matrix< DDRMat > tNodeVal = tGeometryEngine.get_field_vals( tMyGeomIndex,tNewNode.get_index() );
            REQUIRE( tNodeVal(0,0) == Approx( 0.1071 ) );

            //------------------------------------------------------------------------------
            /*
             * ------------------------------------------------------------
             * 4) determine the intersection points along edges [1] and [4]
             * ------------------------------------------------------------
             */

            // edge [1]:
            Matrix< DDRMat > tGlobalPos = {{0},{1}};    // edge [1] goes form x=0 to x=1
            Matrix< DDRMat > tTHat = {{0},{1}};
            Matrix< DDRMat > tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(1)(0,0) }};

            Intersection_Object_Line tIntersectionObject;

            tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
            moris_index tXInd = tGeometryEngine.compute_intersection( tMyGeomIndex, &tIntersectionObject );

            Matrix< F31RMat > tIntersectionAlongX = tGeometryEngine.get_intersection_point( tMyGeomIndex, &tIntersectionObject, tXInd );

            // edge [4]:
            tGlobalPos = {{0},{1}};     // edge [2] goes from y=0 to y=1
            tTHat = {{0},{1}};
            tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(3)(0,0) }};

            tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
            moris_index tYInd = tGeometryEngine.compute_intersection( tMyGeomIndex, &tIntersectionObject );

            Matrix< F31RMat > tIntersectionAlongY = tGeometryEngine.get_intersection_point( tMyGeomIndex, &tIntersectionObject, tYInd );
            //------------------------------------------------------------------------------

//            tIntersectionObject.compute_intersection_sensitivity();
//
//            Matrix< DDRMat > tTemp = tIntersectionObject.get_field_sensitivity_vals( 0 );
            //------------------------------------------------------------------------------
            /*
             * --------------------------------------------------------
             * check determined values
             * --------------------------------------------------------
             */
            CHECK( equal_to( tLSVals(0)(0,0), -0.60 ) );
            CHECK( equal_to( tLSVals(1)(0,0),  0.40 ) );
            CHECK( equal_to( tLSVals(2)(0,0),  (std::sqrt(2)-0.6)) );
            CHECK( equal_to( tLSVals(3)(0,0),  0.40 ) );

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

            CHECK( equal_to( tIntersectionAlongX(0,0), 0.2 ) );
            CHECK( equal_to( tIntersectionAlongY(0,0), 0.2 ) );
            //------------------------------------------------------------------------------
//            std::string tOutputFile = "./analytic_functionalities_2D.exo";
//            tInterpMesh1->create_output_mesh(tOutputFile);
            // clean up
            delete tInterpMesh1;
            delete tIntegMesh1;

        }
//------------------------------------------------------------------------------
TEST_CASE("analytic_geom_setup","[GE],[analytic_geom_setup]")
        {
//            /* ------------------------------------------------------------------------------
//             * Step (1): create the mesh
//             * ------------------------------------------------------------------------------
//             */
//            const std::string tFileName = "generated:10x10x10";
//
//            moris::mtk::Scalar_Field_Info<DDRMat> tInnerSphereField;
//            std::string tSphere1FieldName = "innerSphere";
//            tInnerSphereField.set_field_name(tSphere1FieldName);
//            tInnerSphereField.set_field_entity_rank(EntityRank::NODE);
//
//            moris::mtk::Scalar_Field_Info<DDRMat> tOuterSphereField;
//            std::string tSphere2FieldName = "outerSphere";
//            tOuterSphereField.set_field_name(tSphere2FieldName);
//            tOuterSphereField.set_field_entity_rank(EntityRank::NODE);
//
//            moris::mtk::MtkFieldsInfo tFieldsInfo;
//
//            add_field_for_mesh_input(&tInnerSphereField,tFieldsInfo);
//            add_field_for_mesh_input(&tOuterSphereField,tFieldsInfo);
//
//            mtk::MtkMeshData tMeshData;
//            tMeshData.FieldsInfo = &tFieldsInfo;
//
//            // create mesh pair
//            mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tFileName, &tMeshData );
//            mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh(MeshType::STK, tInterpMesh1);
//
//            // place the pair in mesh manager
//            mtk::Mesh_Manager tMeshManager;
//            uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);
//            /* ------------------------------------------------------------------------------
//             * Step (2): define the geometry representations
//             * ------------------------------------------------------------------------------
//             */
//            moris::Cell< real > tInnerSphereConsts = {{5.0},{5.0},{5.0},{3.0}};
//            moris::Cell< real > tOuterSphereConsts = {{5.0},{5.0},{5.0},{3.5}};
//
//            Ge_Factory tFactory;
//            std::shared_ptr< Geometry > tInnerSphere = tFactory.set_geometry_type(GeomType::ANALYTIC);
//            std::shared_ptr< Geometry > tOuterSphere = tFactory.set_geometry_type(GeomType::ANALYTIC);
//
//            tInnerSphere->set_analytical_function(AnalyticType::SPHERE);
//            tInnerSphere->set_analytical_function_dphi_dx(AnalyticType::SPHERE);
//            tInnerSphere->set_my_mesh(&tMeshManager);
//            tInnerSphere->set_my_constants(tInnerSphereConsts);
//
//            tOuterSphere->set_analytical_function(AnalyticType::SPHERE);
//            tOuterSphere->set_analytical_function_dphi_dx(AnalyticType::SPHERE);
//            tOuterSphere->set_my_mesh(&tMeshManager);
//            tOuterSphere->set_my_constants(tOuterSphereConsts);
//            /* ------------------------------------------------------------------------------
//             * Step (3): build the geometry engine, access data, add geometries to mesh
//             * ------------------------------------------------------------------------------
//             */
//            GE_Core tGeometryEngine;
//            moris_index tInnerSphereIndex = tGeometryEngine.set_geometry( tInnerSphere, tMeshIndex );
//            moris_index tOuterSphereIndex = tGeometryEngine.set_geometry( tOuterSphere, tMeshIndex );
//
//            uint tNumNodes = tMeshManager.get_interpolation_mesh(tMeshIndex)->get_num_nodes();
//            Matrix< DDRMat > tInnerFieldData(tNumNodes,1, 0.0);
//            Matrix< DDRMat > tOuterFieldData(tNumNodes,1, 0.0);
//            for(moris_index n=0; n<(moris_index)tNumNodes; n++)
//            {
//                Matrix< DDRMat > tInnerSphereVals = tGeometryEngine.get_field_vals(tInnerSphereIndex,n);
//                Matrix< DDRMat > tOuterSphereVals = tGeometryEngine.get_field_vals(tOuterSphereIndex,n);
//
//                tInnerFieldData(n) = tInnerSphereVals(0,0);
//                tOuterFieldData(n) = tOuterSphereVals(0,0);
//            }
//            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tSphere1FieldName, EntityRank::NODE, tInnerFieldData);
//            tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tSphere2FieldName, EntityRank::NODE, tOuterFieldData);
//            /* ------------------------------------------------------------------------------
//             * Step (4): build the geometry intersection object and check for an intersection
//             * ------------------------------------------------------------------------------
//             */
//            mtk::Geometry_Type       aGeomType        = mtk::Geometry_Type::LINE;
//            fem::Interpolation_Type  aInterpType      = fem::Interpolation_Type::LAGRANGE;
//            mtk::Interpolation_Order aInterpOrder     = mtk::Interpolation_Order::LINEAR;
//            fem::Interpolation_Type  aTimeInterpType  = fem::Interpolation_Type::LAGRANGE;
//            mtk::Interpolation_Order aTimeInterpOrder = mtk::Interpolation_Order::LINEAR;
//
//            fem::Interpolation_Type  aInterpTypeField  = fem::Interpolation_Type::LAGRANGE;
//            mtk::Interpolation_Order aInterpOrderField = mtk::Interpolation_Order::LINEAR;
//
//            // with the current geometry interpolator available, can only interpolate along one cardinal direction at a time (performing along X1 here)
//            Matrix< DDRMat > tEnd01 = {{ 5,5,5}};
//            Matrix< DDRMat > tEnd02 = {{10,5,5}};
//
//            Matrix<DDRMat> tGlobalPos = {{5},{10}};
//            Matrix< DDRMat > tTHat = {{0},{1}};
//            Matrix< DDRMat > tInnerUHat = {{ tGeometryEngine.compute_value_at_point(tEnd01,tInnerSphereIndex) },{ tGeometryEngine.compute_value_at_point(tEnd02,tInnerSphereIndex) }};
//
//            Intersection_Object tIntObj( aGeomType, aInterpType, aInterpOrder, aTimeInterpType, aTimeInterpOrder, aInterpTypeField, aInterpOrderField );
//            tIntObj.set_coords_and_param_point( tGlobalPos, tTHat, tInnerUHat );
//
//            Matrix< DDRMat > tIntersection = tGeometryEngine.determine_intersection(tInnerSphereIndex, &tIntObj);
//            REQUIRE( tIntersection(0,0) == Approx(8.0));

//            std::string tOutputFile = "./analytic_functionalities_3D.exo";
//            tInterpMesh1->create_output_mesh(tOutputFile);
//            //------------------------------------------------------------------------------
//            delete tInterpMesh1;
//            delete tIntegMesh1;
        }
//------------------------------------------------------------------------------
