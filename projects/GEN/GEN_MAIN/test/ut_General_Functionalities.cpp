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
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_User_Defined_Property.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Phase_Table.hpp"

// HMR includes ---------------------------------
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Field.hpp"

// MTK includes ---------------------------------
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

// XTK include ----------------------------------
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

#include "cl_PRM_HMR_Parameters.hpp"

//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
//-----------------------------------------------------------------------------
real tTempCircleLS( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 0.6;
}
//-----------------------------------------------------------------------------
real tSimpleFunc(const Matrix<DDRMat>& aCoordinates, const Cell<real*>& aPropertyVariables)
{
    // f(x,y) = 2*x + 3*(y*y)
    return (2 * aCoordinates(0) + 3 * aCoordinates(1) * aCoordinates(1));
}

void tSimpleFunc2(const Matrix<DDRMat>& aCoordinates, const Cell<real*>& aPropertyVariables, Matrix<DDRMat>& aSensitivities)
{
}
//-----------------------------------------------------------------------------

//TEST_CASE("general_test_00","[GE],[geom_field_functionality_check]")
//{
//    if(par_size()<=1)
//    {
//        size_t tModelDimension  = 2;
//        uint tLagrangeMeshIndex = 0;
//        //  HMR Parameters setup
//        moris::ParameterList tParameters = prm::create_hmr_parameter_list();
//
//        tParameters.set( "number_of_elements_per_dimension", std::string("2, 2") );
//        tParameters.set( "domain_dimensions",                std::string("2, 2") );
//        tParameters.set( "domain_offset",                    std::string("0, 0") );
//        tParameters.set( "domain_sidesets",                  std::string("1, 2, 3, 4") );
//
//        tParameters.set( "truncate_bsplines", 1 );
//        tParameters.set( "lagrange_orders",   std::string("1") );
//        tParameters.set( "lagrange_pattern",  std::string("0") );
//        tParameters.set( "bspline_orders",    std::string("1") );
//        tParameters.set( "bspline_pattern",   std::string("0") );
//
//        tParameters.set( "lagrange_output_meshes", std::string("0") );
//        tParameters.set( "lagrange_input_meshes",  std::string("0") );
//
//        tParameters.set( "lagrange_to_bspline", std::string("0") );
//
//        tParameters.set( "use_multigrid", 0 );
//
//        tParameters.set( "refinement_buffer", 1 );
//        tParameters.set( "staircase_buffer",  1 );
//
//        tParameters.insert( "initial_refinement", 0 );
//
//        //  HMR Initialization
//        moris::hmr::HMR tHMR( tParameters );
//
//        auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >
//
//        tHMR.perform_initial_refinement( 0 );
//        //------------------------------------------------------------------------------
//        tDatabase->update_bspline_meshes();
//        tDatabase->update_lagrange_meshes();
//        //------------------------------------------------------------------------------
//        tHMR.finalize();
//
//        hmr::Interpolation_Mesh_HMR *      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//        moris::hmr::Integration_Mesh_HMR * tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );
//
//        mtk::Mesh_Manager tMesh;
//        moris_index tMeshPairIndex = tMesh.register_mesh_pair( tInterpMesh, tIntegrationMesh );
//        //------------------------------------------------------------------------------
//        moris::Cell<moris::mtk::Vertex const *> tAllVertices = tMesh.get_interpolation_mesh( tMeshPairIndex )->get_all_vertices();
//
//        uint tNumVertices = tAllVertices.size();
//
//        moris::Matrix< DDRMat > tCoords(tNumVertices,tModelDimension);
//
//        for(uint i=0; i<tNumVertices; i++)  // build a matrix with all node coordinates
//        {
//            tCoords.set_row( i, tAllVertices(i)->get_coords() );
//        }
//
//        Matrix<DDRMat> tAdvs(1, 1, 0.0);
//        std::shared_ptr< Property > tFieldProperty = std::make_shared< User_Defined_Property >(tAdvs,
//                                                                                               Matrix<DDUMat>(0, 0),
//                                                                                               Matrix<DDUMat>(0, 0),
//                                                                                               Matrix<DDRMat>(0, 0),
//                                                                                               Cell<std::shared_ptr<Property>>(0),
//                                                                                               &tSimpleFunc,
//                                                                                               &tSimpleFunc2);
//
//        tField.initialize( &tMesh );
//
//        Cell<std::shared_ptr<ge::Geometry_Discrete>> tGeometryVector(1);
//        tGeometryVector(0) = std::make_shared<moris::ge::Geometry_Field>(&tField);
//
//        moris::ge::Phase_Table      tPhaseTable( tGeometryVector.size(), moris::ge::Phase_Table_Structure::EXP_BASE_2 );
//        moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tModelDimension );
//        //=================== manual calls to GE (w/out XTK model) =============================
//        tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumVertices );
//
//        Matrix< DDRMat > tLSVals( tNumVertices,1 );
//        for( uint i=0; i<tNumVertices; i++ )
//        {
//            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
//            tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );
//        }
//        //============================= perform checks =========================================
//        real tEpsilon = 1e-8;
//        Matrix< DDRMat > tCorrectVals = { {  0.0 },
//                {  2.0 },
//                {  5.0 },
//                {  3.0 },
//                {  4.0 },
//                {  7.0 },
//                { 14.0 },
//                { 12.0 },
//                { 16.0 } };
//
//        for( uint i=0; i<tNumVertices; i++ )
//        {
//            REQUIRE( std::abs(tLSVals(i) - tCorrectVals(i)) <tEpsilon );
//        }
//        delete tInterpMesh;
//        delete tIntegrationMesh;
//        //================================ end =================================================
//    }   // end serial statement
//}

//-----------------------------------------------------------------------------
TEST_CASE("general_test_01","[GE],[sensitivity_check_01]")
{
    if(par_size()<=1)
    {
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
        mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
        mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
        //------------------------------------------------------------------------------
        real tRadius = 0.6;
        Cell<std::shared_ptr<moris::ge::Geometry_Analytic>> tGeometry(1);
        tGeometry(0) = std::make_shared<moris::ge::Circle>(0.0, 0.0, tRadius);

        moris::ge::Phase_Table      tPhaseTable( 1, moris::ge::Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometry, tPhaseTable, tNumDim );

        //=================== manual calls to GE (w/out XTK model) =============================
        uint tNumNodes = tInterpMesh->get_num_nodes();
        tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );

        Matrix< DDRMat > tLSVals( tNumNodes,1 );
        for( uint i=0; i<tNumNodes; i++ )
        {
            tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
            tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );
        }
        tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);

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
        //======================================
        REQUIRE( tNumIntersections == 2 );  // there should be two intersected edges
        //======================================
        Matrix< DDRMat >   tLclCoord( tNumIntersections,1 );
        Matrix< DDRMat >   tNewNodeCoords( tNumIntersections,tNumDim );
        for( uint i=0; i<tNumIntersections; i++ )
        {
            tLclCoord(i) = tGeometryObjects(i).get_interface_lcl_coord();
            tNewNodeCoords.set_row( i,tGeometryObjects(i).get_interface_glb_coord() );
        }
        //================== perform checks before moving on ====================
        REQUIRE( tLclCoord(0) =  0.2 );
        REQUIRE( tLclCoord(1) = -0.2 );
        REQUIRE( tNewNodeCoords(0,0) = 0.6 );
        REQUIRE( tNewNodeCoords(1,1) = 0.6 );
        //=======================================================================
        Matrix< moris::IndexMat > tEdgeIndices00{{ 0,1 }};
        std::shared_ptr< xtk::Topology > tTop00 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices00 );

        Matrix< moris::IndexMat > tEdgeIndices01{{ 0,3 }};
        std::shared_ptr< xtk::Topology > tTop01 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices01 );
        moris::Cell< xtk::Topology* > tTopCell(2);
        tTopCell(0) = tTop00.get();          tTopCell(1) = tTop01.get();

        moris::Cell< moris_index > tIndices(2);
        tIndices(0) = 4;   tIndices(1) = 5;
        moris::Cell< Matrix<DDRMat> > tLocalCoords(2);
        tLocalCoords(0) = {{0.2}};  tLocalCoords(1) = {{-0.2}};
        moris::Cell< Matrix<DDRMat> > tGlbCoords(2);
        tGlbCoords(0) = tNewNodeCoords.get_row(0); tGlbCoords(1) = tNewNodeCoords.get_row(1);

        tGENGeometryEngine.create_new_node_geometry_objects( tIndices, true, tTopCell, tLocalCoords, tGlbCoords );

        Matrix< DDRMat > tAllCoords = {    { 0.0, 0.0 },
                                           { 1.0, 0.0 },
                                           { 1.0, 1.0 },
                                           { 0.0, 1.0 },
                                       { tRadius, 0.0 },
                                       { 0.0, tRadius }};

        Matrix< IndexMat > tInd{{ tIndices(0),tIndices(1) }};
        tGENGeometryEngine.compute_interface_sensitivity( tInd, tAllCoords, (moris_index)0 );

        GEN_Geometry_Object tGeomObj01 = tGENGeometryEngine.get_geometry_object( tIndices(0) );
        moris::Matrix< moris::DDRMat >tSens01 = tGeomObj01.get_sensitivity_dx_dp();
//print( tSens01,"tSens01" );

        GEN_Geometry_Object tGeomObj02 = tGENGeometryEngine.get_geometry_object( tIndices(1) );
        moris::Matrix< moris::DDRMat >tSens02 = tGeomObj02.get_sensitivity_dx_dp();
//print( tSens02,"tSens02" );
        //================================== end ===============================================

        //------------------------------------------------------------------------------
//        std::string tOutputFile = "./output_sensitivityCheck_01.exo";
//        tInterpMesh->create_output_mesh(tOutputFile);
        delete tInterpMesh;
        delete tIntegMesh;
        //------------------------------------------------------------------------------
    }
}

//-----------------------------------------------------------------------------
TEST_CASE("general_test_02","[GE],[sensitivity_check_02]")
{
    if(par_size()<=1)
    {
        size_t tModelDimension = 2;
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        moris::ParameterList tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", std::string("1, 1") );
        tParameters.set( "domain_dimensions",                std::string("1, 1") );
        tParameters.set( "domain_offset",                    std::string("0, 0") );

//        tParameters.set( "domain_sidesets",      "1, 2, 3, 4, 5, 6" );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_output_meshes", std::string("0") );
        tParameters.set( "lagrange_input_meshes", std::string("0") );

        tParameters.set( "lagrange_to_bspline", std::string("0") );

        tParameters.set( "use_multigrid", 0 );

        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );

        tParameters.insert( "initial_refinement", 0 );

        //  HMR Initialization
        moris::hmr::HMR tHMR( tParameters );

        std::shared_ptr< hmr::Database > tDatabase = tHMR.get_database();

        tHMR.perform_initial_refinement( 0 );

//        tDatabase->update_bspline_meshes();
//        tDatabase->update_lagrange_meshes();

        std::shared_ptr< moris::hmr::Mesh > tHMRMesh = tHMR.create_mesh( tLagrangeMeshIndex );
//        tHMR.finalize();

        hmr::Interpolation_Mesh_HMR * tInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
//        std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );
        //------------------------------------------------------------------------------
        std::shared_ptr< hmr::Field > tField = tHMRMesh->create_field( "circle", tLagrangeMeshIndex);

        tField->evaluate_scalar_function( moris::ge::tTempCircleLS );

        tDatabase->update_bspline_meshes();
        tDatabase->update_lagrange_meshes();

        tHMR.finalize();
//        tHMR.save_to_exodus( 0, "output_sensitivityCheck_02.g" );
        //------------------------------------------------------------------------------
        Cell<std::shared_ptr<moris::ge::Geometry_Analytic>> tGeometry(1);
        tGeometry(0) = std::make_shared<moris::ge::Circle>(0.0, 0.0, 0.6);

        moris::ge::Phase_Table      tPhaseTable( 1,  moris::ge::Phase_Table_Structure::EXP_BASE_2 );
        moris::ge::Geometry_Engine  tGENGeometryEngine( tGeometry, tPhaseTable, tModelDimension );

        //------------------------------------------------------------------------------
        xtk::Model tXTKModel( tModelDimension, tInterpMesh, &tGENGeometryEngine );
        tXTKModel.mVerbose = false;
        Cell< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4 };
        tXTKModel.decompose( tDecompositionMethods );

        tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

        xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh( );
        xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh( );

        // place the pair in mesh manager
        mtk::Mesh_Manager tMeshManager;
        tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh);

        xtk::Output_Options tOutputOptions;
        tOutputOptions.mAddNodeSets = false;
        tOutputOptions.mAddSideSets = false;
        tOutputOptions.mAddClusters = false;

        std::string tCircleFieldName = "circle";
        tOutputOptions.mRealNodeExternalFieldNames = { tCircleFieldName };

        moris::mtk::Integration_Mesh* tIntegMesh1 = tXTKModel.get_output_mesh( tOutputOptions );
        //------------------------------------------------------------------------------
//        Matrix< DDRMat > tTMatrix = tHMRMesh->get_t_matrix_of_node_loc_ind( 0, EntityRank::ELEMENT );
//        print( tTMatrix,"tTMatrix" );

        //------------------------------------------------------------------------------
        delete tIntegMesh1;
        delete tInterpMesh;

    }
}

//-----------------------------------------------------------------------------

}   // ge namespace
}   // moris namespace
