/*
 * cl_Pdv_Tests.cpp
 *
 *  Created on: Dec 20, 2019
 *      Author: sonne
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

// GE include -----------------------------------
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Field_User_Defined.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Phase_Table.hpp"
#include "cl_GEN_Property.hpp"

// HMR includes ---------------------------------
#include "cl_HMR.hpp"
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
//------------------------------------------------------------------------------

namespace moris
{
    namespace ge
    {
    Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff,
                                        moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator )
    {
        return aCoeff( 0 );
    }
    //------------------------------------------------------------------------------

    TEST_CASE("pdv_test_00","[GE],[pdv_check_00]")
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
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField00;
            std::string tField00Name = "circle";
            tNodeField00.set_field_name(tField00Name);
            tNodeField00.set_field_entity_rank(EntityRank::NODE);

            // declare scalar node field for the density
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField01;
            std::string tField01Name = "density";
            tNodeField01.set_field_name(tField01Name);
            tNodeField01.set_field_entity_rank(EntityRank::NODE);

            // declare scalar node field for the temperature
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField02;
            std::string tField02Name = "temperature";
            tNodeField02.set_field_name(tField02Name);
            tNodeField02.set_field_entity_rank(EntityRank::NODE);

            // initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;
            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField00,tFieldsInfo);
            add_field_for_mesh_input(&tNodeField01,tFieldsInfo);
            add_field_for_mesh_input(&tNodeField02,tFieldsInfo);

            // declare some supplementary fields
            tMeshData.FieldsInfo = &tFieldsInfo;

            // create mesh pair
            mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
            mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

            // place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
            //------------------------------------------------------------------------------
            //------------------------------------------------------------------------------
            Cell< enum GEN_PDV >  tPdvList(2);
            tPdvList(0) = GEN_PDV::DENSITY;
            tPdvList(1) = GEN_PDV::TEMP;

            std::shared_ptr< GEN_Property > tConstDensityProp = std::make_shared< GEN_Property >();
            tConstDensityProp->set_parameters( { {{ 1234 }} } );
            tConstDensityProp->set_val_function( tConstValFunction );

            std::shared_ptr< GEN_Property > tConstTempProp = std::make_shared< GEN_Property >();
            tConstTempProp->set_parameters( { {{ 99 }} } );
            tConstTempProp->set_val_function( tConstValFunction );

            Cell< GEN_Property* > tPropertyList(2);
            tPropertyList(0) = tConstDensityProp.get();
            tPropertyList(1) = tConstTempProp.get();

            //------------------------------------------------------------------------------
            real tRadius  = 0.6;
            real tXcenter = 0;
            real tYcenter = 0;
            moris::ge::Circle tCircle( tRadius, tXcenter, tYcenter );
            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tCircle };

            moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(), Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tNumDim );

            uint tNumNodes = tInterpMesh->get_num_nodes();
            tGENGeometryEngine.set_pdv_property_list( tPdvList, tPropertyList );
            //=================== manual calls to GE (w/out XTK model) =============================

            tGENGeometryEngine.initialize_pdv_hosts_for_background_mesh_nodes( tNumNodes );
            tGENGeometryEngine.initialize_geometry_objects_for_background_mesh_nodes( tNumNodes );

            Pdv_Host_Manager* tAllPdvHosts = tGENGeometryEngine.get_pdv_hosts();

            Cell< Matrix< DDRMat > > tDensPdvVals( tNumNodes );
            Cell< Matrix< DDRMat > > tTempPdvVals( tNumNodes );

            Matrix< DDRMat > tLSVals( tNumNodes, 1 );
            Matrix< DDRMat > tDensityVals( tNumNodes, 1 );
            Matrix< DDRMat > tTempVals( tNumNodes, 1 );

            for( uint i=0; i<tNumNodes; i++ )
            {
                tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
                tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );

                tAllPdvHosts->get_pdv_values( i, GEN_PDV::DENSITY, tDensPdvVals(i) );
                tDensityVals(i) = tDensPdvVals(i)(0,0);

                tAllPdvHosts->get_pdv_values( i, GEN_PDV::TEMP, tTempPdvVals(i) );
                tTempVals(i) = tTempPdvVals(i)(0,0);
            }

            //----- check node values -----
            for(uint i=0; i<tNumNodes; i++)
            {
                REQUIRE( tDensPdvVals(i)(0,0) = 1234 );
                REQUIRE( tTempPdvVals(i)(0,0) = 99 );
            }
            //-----------------------------
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tField00Name, EntityRank::NODE, tLSVals);         // add LS values to mesh
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tField01Name, EntityRank::NODE, tDensityVals);    // add density values to mesh
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tField02Name, EntityRank::NODE, tTempVals);       // add temperature values to mesh

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
            tTopCell(0) = tTop00.get();    tTopCell(1) = tTop01.get();

            moris::Cell< moris_index > tIndices(2);
            tIndices(0) = 4;    tIndices(1) = 5;
            moris::Cell< Matrix<DDRMat> > tLocalCoords(2);
            tLocalCoords(0) = {{0.2}};    tLocalCoords(1) = {{-0.2}};
            moris::Cell< Matrix<DDRMat> > tGlbCoords(2);
            tGlbCoords(0) = tNewNodeCoords.get_row(0);     tGlbCoords(1) = tNewNodeCoords.get_row(1);

            tGENGeometryEngine.create_new_node_geometry_objects( tIndices, true, tTopCell, tLocalCoords, tGlbCoords );
            //============================ end manual calls ========================================

            //------------------------------------------------------------------------------
//            std::string tOutputFile = "./pdvCheck.exo";
//            tInterpMesh->create_output_mesh(tOutputFile);
            delete tInterpMesh;
            delete tIntegMesh;
            //------------------------------------------------------------------------------
        }
    }

    //------------------------------------------------------------------------------
    TEST_CASE("pdv_test_01","[GE],[pdv_check_01]")
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
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField00;
            std::string tField00Name = "circle";
            tNodeField00.set_field_name(tField00Name);
            tNodeField00.set_field_entity_rank(EntityRank::NODE);

            // declare scalar node field for the density
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField01;
            std::string tField01Name = "density";
            tNodeField01.set_field_name(tField01Name);
            tNodeField01.set_field_entity_rank(EntityRank::NODE);

            // declare scalar node field for the temperature
            moris::mtk::Scalar_Field_Info<DDRMat> tNodeField02;
            std::string tField02Name = "temperature";
            tNodeField02.set_field_name(tField02Name);
            tNodeField02.set_field_entity_rank(EntityRank::NODE);

            // initialize field information container
            moris::mtk::MtkFieldsInfo tFieldsInfo;
            // Place the node field into the field info container
            add_field_for_mesh_input(&tNodeField00,tFieldsInfo);
            add_field_for_mesh_input(&tNodeField01,tFieldsInfo);
            add_field_for_mesh_input(&tNodeField02,tFieldsInfo);

            // declare some supplementary fields
            tMeshData.FieldsInfo = &tFieldsInfo;

            // create mesh pair
            mtk::Interpolation_Mesh* tInterpMesh = create_interpolation_mesh( MeshType::STK, tMeshData );
            mtk::Integration_Mesh*   tIntegMesh  = create_integration_mesh_from_interpolation_mesh( MeshType::STK,tInterpMesh );

            // place the pair in mesh manager
            mtk::Mesh_Manager tMeshManager;
            uint tMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );
            //------------------------------------------------------------------------------
            //------------------------------------------------------------------------------
            Cell< enum GEN_PDV >  tPdvList(2);
            tPdvList(0) = GEN_PDV::DENSITY;
            tPdvList(1) = GEN_PDV::TEMP;

            std::shared_ptr< GEN_Property > tConstDensityProp = std::make_shared< GEN_Property >();
            tConstDensityProp->set_parameters( { {{ 1234 }} } );
            tConstDensityProp->set_val_function( tConstValFunction );

            std::shared_ptr< GEN_Property > tConstTempProp = std::make_shared< GEN_Property >();
            tConstTempProp->set_parameters( { {{ 99 }} } );
            tConstTempProp->set_val_function( tConstValFunction );

            Cell< GEN_Property* > tPropertyList(2);
            tPropertyList(0) = tConstDensityProp.get();
            tPropertyList(1) = tConstTempProp.get();

            //------------------------------------------------------------------------------
            real tRadius  = 0.6;
            real tXcenter = 0;
            real tYcenter = 0;
            moris::ge::Circle tCircle( tRadius, tXcenter, tYcenter );
            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tCircle };

            moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(), Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGENGeometryEngine( tGeometryVector, tPhaseTable, tNumDim );

            tGENGeometryEngine.set_pdv_property_list( tPdvList, tPropertyList );

            xtk::Model tXTKModel( tNumDim, tInterpMesh, tGENGeometryEngine );
            tXTKModel.mVerbose = false;

            Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4};
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( EntityRank::NODE, 0 );

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh( );
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh( );
            uint tEnrMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegMesh );

            // --------- create the pdv objects for the new enriched mesh ---------
            moris::ge::GEN_Geometry_Engine tXTKGeomEng = tXTKModel.get_geom_engine();

            uint tNumNodes = tEnrInterpMesh.get_num_nodes();

            tXTKGeomEng.initialize_pdv_hosts_for_background_mesh_nodes( tNumNodes );

            Pdv_Host_Manager* tAllPdvHosts = tXTKGeomEng.get_pdv_hosts();
            // --------------------------------------------------------------------
            Matrix< IndexMat > tIndices;
            tAllPdvHosts->get_all_node_indices(tIndices);

            Cell< Matrix< DDRMat > > tDensPdvVals( tNumNodes );
            Cell< Matrix< DDRMat > > tTempPdvVals( tNumNodes );

            uint tIters = tIndices.n_rows();
            for( uint i=0; i<tIters; i++ )
            {
                tAllPdvHosts->get_pdv_values( tIndices(i), GEN_PDV::TEMP, tTempPdvVals(i) );
                tAllPdvHosts->get_pdv_values( tIndices(i), GEN_PDV::DENSITY, tDensPdvVals(i) );
            }

            //----- check node values -----
            for(uint i=0; i<tIters; i++)
            {
                REQUIRE( tDensPdvVals(i)(0,0) = 1234 );
                REQUIRE( tTempPdvVals(i)(0,0) = 99 );
            }
            //-----------------------------

            //------------------------------------------------------------------------------
            delete tInterpMesh;
            delete tIntegMesh;
        }
    }

    //------------------------------------------------------------------------------

    }   // end ge namespace
}       //end moris namespace



