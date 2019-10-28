/*
 * UT_Finite_Difference_Checks.cpp
 *
 *  Created on: Sep 10, 2019
 *      Author: sonne
 */
#include <catch.hpp>

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Intersection_Object_Line.hpp"

// MTK includes
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

// LINALG includes
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

using namespace moris;
using namespace mtk;
using namespace ge;

TEST_CASE("Finite_Difference_Checks_01","[moris],[GE],[FD_on_circle]")
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

        // input parameters for the circle LS
        moris::Cell< real > tCircleInputs = {{0.6},
                                             {0},
                                             {0}};
        //------------------------------------------------------------------------------

        Ge_Factory tFactory;
        std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(GeomType::ANALYTIC);

        tGeom1->set_my_mesh(&tMeshManager);

        moris_index tSubIndex = tGeom1->set_analytical_function( AnalyticType::CIRCLE, tCircleInputs );

        GE_Core tGeometryEngine;
        moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );

        PDV_Info* tPDVInfo = tGeometryEngine.get_pdv_info_pointer( tMyGeomIndex );

        Matrix< DDRMat > tLSVals                       = tPDVInfo->get_field_vals( tSubIndex );           // phi
        moris::Cell< Matrix< DDRMat > > tSensitivities = tPDVInfo->get_sensitivity_vals( tSubIndex );     // dphi/dp

//        tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);      // add the determined values as a field on the mesh (for output purposes)

        // edge [1]:
        Matrix< DDRMat > tGlobalPos = {{0,0},
                                       {1,0}};
        Matrix< DDRMat > tTHat = {{0},
                                  {1}};
        Matrix< DDRMat > tUHat = {{ tLSVals(0) },
                                  { tLSVals(1) }};

        Intersection_Object_Line tIntersectionObject;

        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
        /*
         * (1) compute and store xgamma, dxgamma/dphi
         */
        moris_index tXInd = tPDVInfo->compute_intersection( &tIntersectionObject );

        Matrix< F31RMat > tXGamma = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd );

        tPDVInfo->compute_intersection_sensitivity( &tIntersectionObject, tXInd );

        Matrix< DDRMat > tdXGamma_dphi = tPDVInfo->get_intersection_sensitivity( &tIntersectionObject, tXInd );
        //------------------------------------------------------------------------------
        /*
         * perform perturbations on the first field value (phi_A) and compute finite difference
         */
        moris::real epsilon = 0.000001;
        /*
         * (2) perturb phi up
         */
        tUHat = {{ tLSVals(0)+epsilon },{ tLSVals(1) }};
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
        moris_index tXInd_up = tPDVInfo->compute_intersection( &tIntersectionObject );
        Matrix< F31RMat > tGlobPoint_up = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_up );
        /*
         * (3) perturb phi down
         */
        tUHat = {{ tLSVals(0)-epsilon },{ tLSVals(1) }};
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
        moris_index tXInd_down = tPDVInfo->compute_intersection( &tIntersectionObject );
        Matrix< F31RMat > tGlobPoint_down = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_down );
        //------------------------------------------------------------------------------
        /*
         * (4) compute dxgamma/dphi (FD)
         */
        Matrix< F31RMat > tdXGamma_dphiA_FD = (tGlobPoint_up-tGlobPoint_down)/(2*epsilon);
        /*
         * --------------------------------------------------------
         * check determined values
         * --------------------------------------------------------
         */
        REQUIRE( tdXGamma_dphiA_FD(0,0) == Approx(tdXGamma_dphi(0,0)));
        REQUIRE( tdXGamma_dphiA_FD(0,1) == Approx(tdXGamma_dphi(0,1)));

        //------------------------------------------------------------------------------
        /*
         * perform perturbations on the second field value (phi_B) and compute finite difference
         */
        tUHat = {{ tLSVals(0) },{ tLSVals(1)+epsilon }};
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
        tXInd_up      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_up = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_up );

        tUHat = {{ tLSVals(0) },{ tLSVals(1)-epsilon }};
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
        tXInd_down      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_down = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_down );

        Matrix< F31RMat > tdXGamma_dphiB_FD = (tGlobPoint_up-tGlobPoint_down)/(2*epsilon);

        REQUIRE( tdXGamma_dphiB_FD(0,0) == Approx(tdXGamma_dphi(1,0)));
        REQUIRE( tdXGamma_dphiB_FD(0,1) == Approx(tdXGamma_dphi(1,1)));
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        /*
         * FD check on the pdv sensitivity dxgamma/dp
         */
        Matrix< DDRMat > tdXGamma_dp = tPDVInfo->get_dxgamma_dp( &tIntersectionObject );    // initial dxgamma/dp
        /*
         * purturb the first pdv up and recalculate xgamma
         */
        tCircleInputs = {{0.6+epsilon},
                         {0},
                         {0}};
        tGeom1->set_constants( tCircleInputs, tSubIndex );
        tPDVInfo->reinitialize_data_table( tGeom1, tSubIndex );

        tLSVals = tPDVInfo->get_field_vals( tSubIndex );

        tUHat = {{ tLSVals(0) },{ tLSVals(1) }};    // correct field values
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );

        tXInd_up      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_up = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_up );
        /*
         * purturb the first pdv down and recalculate xgamma
         */
        tCircleInputs = {{0.6-epsilon},
                         {0},
                         {0}};
        tGeom1->set_constants( tCircleInputs, tSubIndex );
        tPDVInfo->reinitialize_data_table( tGeom1, tSubIndex );

        tLSVals = tPDVInfo->get_field_vals( tSubIndex );

        tUHat = {{ tLSVals(0) },{ tLSVals(1) }};    // correct field values
        tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );

        tXInd_down      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_down = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_down );
        //------------------------------------------------------------------------------
        Matrix< F31RMat > tdXGamma_dp_FD = (tGlobPoint_up-tGlobPoint_down)/(2*epsilon);

        REQUIRE( tdXGamma_dp_FD(0,0) == Approx(-tdXGamma_dp(0,0)) );
        //------------------------------------------------------------------------------
        //TODO add test for y-intersection point
//        std::string tOutputFile = "./ge_finite_diff_ut.exo";
//        tInterpMesh1->create_output_mesh(tOutputFile);
        // clean up
        delete tInterpMesh1;
        delete tIntegMesh1;
}
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

TEST_CASE("Finite_Difference_Checks_02","[moris],[GE],[FD_on_sphere]")
{
    if(par_size()<=1)
    {
        uint aNumElemTypes = 1;     // hex
        uint aNumDim = 3;           // specify number of spatial dimensions

        Matrix< IdMat > aElementConnHex = {{ 1, 2, 3, 4, 5, 6, 7, 8 }};   // specify element connectivity of hex for mesh

        Matrix< IdMat > aElemLocalToGlobalHex = {{ 1 }};      // specify the local to global element map for hexs

        Matrix< DDRMat > aCoords = {{ 0.0, 0.0, 0.0 },
                                    { 1.0, 0.0, 0.0 },
                                    { 1.0, 0.0, 1.0 },
                                    { 0.0, 0.0, 1.0 },

                                    { 0.0, 1.0, 0.0 },
                                    { 1.0, 1.0, 0.0 },
                                    { 1.0, 1.0, 1.0 },
                                    { 0.0, 1.0, 1.0 }};             // Node coordinate matrix

        Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4, 5, 6, 7, 8 }}; // specify the local to global map
        //------------------------------------------------------------------------------
        // create MORIS mesh using MTK database
        mtk::MtkMeshData tMeshData( aNumElemTypes );
        tMeshData.CreateAllEdgesAndFaces = true;
        tMeshData.SpatialDim = & aNumDim;
        tMeshData.ElemConn(0) = & aElementConnHex;
        //------------------------------------------------------------------------------
        tMeshData.NodeCoords = & aCoords;
        tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalHex;
        tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
        //------------------------------------------------------------------------------
        // declare scalar node field for the circle LS
        moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
        std::string tFieldName = "sphere";
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

        // input parameters for the circle LS
        moris::Cell< real > tSphereInputs = {{0.6},
                                             {0.0},
                                             {0.0},
                                             {0.6}};
        //------------------------------------------------------------------------------

        Ge_Factory tFactory;
        std::shared_ptr< Geometry > tGeom = tFactory.set_geometry_type(GeomType::ANALYTIC);

        tGeom->set_my_mesh( &tMeshManager );

        moris_index tSubIndex = tGeom->set_analytical_function( AnalyticType::SPHERE, tSphereInputs );

        GE_Core tGeometryEngine;
        moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom );

        PDV_Info* tPDVInfo = tGeometryEngine.get_pdv_info_pointer( tMyGeomIndex );

        Matrix< DDRMat > tLSVals                       = tPDVInfo->get_field_vals( tSubIndex );           // phi
        moris::Cell< Matrix< DDRMat > > tSensitivities = tPDVInfo->get_sensitivity_vals( tSubIndex );     // dphi/dp

        tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);      // add the determined values as a field on the mesh (for output purposes)

        // edge [1]:
        Matrix< DDRMat > tGlobalPos(2,3);
        tGlobalPos.get_row( 0 ) = tGeom->get_my_mesh()->get_integration_mesh( 0 )->get_node_coordinate( 0 ).get_row(0);
        tGlobalPos.get_row( 1 ) = tGeom->get_my_mesh()->get_integration_mesh( 0 )->get_node_coordinate( 1 ).get_row(0);

        Matrix< DDRMat > tTHat = {{0},
                                  {1}};
        Matrix< DDRMat > tUHat = {{ tLSVals(0) },
                                  { tLSVals(1) }};

        Intersection_Object_Line tIntersectionObject;

        tIntersectionObject.set_coords_and_param_point( tGeom, tGlobalPos, tTHat, tUHat );
        /*
         * (1) compute and store xgamma, dxgamma/dphi
         */
        moris_index tXInd = tPDVInfo->compute_intersection( &tIntersectionObject );

        Matrix< F31RMat > tXGamma = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd );

        tPDVInfo->compute_intersection_sensitivity( &tIntersectionObject, tXInd );

        Matrix< DDRMat > tdXGamma_dphi = tPDVInfo->get_intersection_sensitivity( &tIntersectionObject, tXInd );
        //------------------------------------------------------------------------------
        /*
         * perform perturbations on the first field value and compute finite difference
         */
        moris::real epsilon = 0.000001;
        /*
         * (2) perturb phi up
         */
        tUHat = {{ tLSVals(0)+epsilon },{ tLSVals(1) }};
        tIntersectionObject.set_coords_and_param_point( tGeom, tGlobalPos, tTHat, tUHat );
        moris_index tXInd_up = tPDVInfo->compute_intersection( &tIntersectionObject );
        Matrix< F31RMat > tGlobPoint_up = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_up );
        /*
         * (3) perturb phi down
         */
        tUHat = {{ tLSVals(0)-epsilon },{ tLSVals(1) }};
        tIntersectionObject.set_coords_and_param_point( tGeom, tGlobalPos, tTHat, tUHat );
        moris_index tXInd_down = tPDVInfo->compute_intersection( &tIntersectionObject );
        Matrix< F31RMat > tGlobPoint_down = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_down );
        //------------------------------------------------------------------------------
        /*
         * (4) compute dxgamma/dphi (FD)
         */
        Matrix< F31RMat > tdXGamma_dphi_FD = (tGlobPoint_up-tGlobPoint_down)/(2*epsilon);
        /*
         * --------------------------------------------------------
         * check determined values
         * --------------------------------------------------------
         */
        REQUIRE( tdXGamma_dphi_FD(0,0) == Approx(tdXGamma_dphi(0,0)));
        REQUIRE( tdXGamma_dphi_FD(0,1) == Approx(tdXGamma_dphi(0,1)));
        //------------------------------------------------------------------------------
        /*
         * perform perturbations on the second field value and compute finite difference
         */
        tUHat = {{ tLSVals(0) },
                 { tLSVals(1)+epsilon }};
        tIntersectionObject.set_coords_and_param_point( tGeom, tGlobalPos, tTHat, tUHat );
        tXInd_up      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_up = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_up );

        tUHat = {{ tLSVals(0) },
                 { tLSVals(1)-epsilon }};
        tIntersectionObject.set_coords_and_param_point( tGeom, tGlobalPos, tTHat, tUHat );
        tXInd_down      = tPDVInfo->compute_intersection( &tIntersectionObject );
        tGlobPoint_down = tPDVInfo->get_intersection_point_global_coord( &tIntersectionObject, tXInd_down );

        tdXGamma_dphi_FD = (tGlobPoint_up-tGlobPoint_down)/(2*epsilon);

        REQUIRE( tdXGamma_dphi_FD(0,0) == Approx(tdXGamma_dphi(1,0)));
        REQUIRE( tdXGamma_dphi_FD(0,1) == Approx(tdXGamma_dphi(1,1)));
        //------------------------------------------------------------------------------

//        std::string tOutputFile = "./ge_finite_diff_ut.exo";
//        tInterpMesh1->create_output_mesh(tOutputFile);
        // clean up
        delete tInterpMesh1;
        delete tIntegMesh1;
    }
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------


