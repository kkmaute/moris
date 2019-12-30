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
                                        moris::Cell< fem::Field_Interpolator* > & aDvFieldInterpolator,
                                        fem::Geometry_Interpolator              * aGeometryInterpolator )
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
            //==============================================================================

            Cell< enum GEN_PDV >  tPdvList{{ GEN_PDV::DENSITY }};

            std::shared_ptr< GEN_Property > tConstDensityProp = std::make_shared< GEN_Property >();
            tConstDensityProp->set_parameters( { {{ 98.2 }} } );
            tConstDensityProp->set_val_function( tConstValFunction );

            Cell< GEN_Property* > tPropertyList(1);
            tPropertyList(0) = tConstDensityProp.get();

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

            tGENGeometryEngine.set_pdv_property_list( tPdvList, tPropertyList );
            tGENGeometryEngine.initialize_pdv_hosts_for_background_mesh_nodes( tNumNodes );

            Pdv_Host_Manager tAllPdvHosts = tGENGeometryEngine.get_pdv_hosts();

            Cell< Matrix< DDRMat > > tPdvVals( tNumNodes );
            Matrix< DDRMat > tLSVals( tNumNodes,1 );
            for( uint i=0; i<tNumNodes; i++ )
            {
                tGENGeometryEngine.initialize_geometry_object_phase_values( tCoords );
                tLSVals(i) = tGENGeometryEngine.get_entity_phase_val( i,0 );

                tAllPdvHosts.get_pdv_host_from_manager(i).get_pdv_values( GEN_PDV::DENSITY, tPdvVals(i) );
print( tPdvVals(i),"tPdvVals(i)" );
            }
            tInterpMesh->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);



//            moris::Matrix< moris::IndexMat > tNodetoEdgeConnectivity(4, 2);
//            // Edge 0
//            (tNodetoEdgeConnectivity)(0, 0) = 0;
//            (tNodetoEdgeConnectivity)(0, 1) = 1;
//            // Edge 1
//            (tNodetoEdgeConnectivity)(1, 0) = 1;
//            (tNodetoEdgeConnectivity)(1, 1) = 2;
//            // Edge 2
//            (tNodetoEdgeConnectivity)(2, 0) = 2;
//            (tNodetoEdgeConnectivity)(2, 1) = 3;
//            // Edge 3
//            (tNodetoEdgeConnectivity)(3, 0) = 3;
//            (tNodetoEdgeConnectivity)(3, 1) = 0;
//
//            Cell< GEN_Geometry_Object > tGeometryObjects;
//            tGENGeometryEngine.is_intersected( tCoords, tNodetoEdgeConnectivity, (size_t) 1, tGeometryObjects );
//
//            size_t tNumIntersections = tGeometryObjects.size();
//            //======================================
//            REQUIRE( tNumIntersections == 2 );  // there should be two intersected edges
//            //======================================
//            Matrix< DDRMat >   tLclCoord( tNumIntersections,1 );
//            Matrix< DDRMat >   tNewNodeCoords( tNumIntersections,tNumDim );
//            for( uint i=0; i<tNumIntersections; i++ )
//            {
//                tLclCoord(i) = tGeometryObjects(i).get_interface_lcl_coord();
//                tNewNodeCoords.set_row( i,tGeometryObjects(i).get_interface_glb_coord() );
//            }
//            //================== perform checks before moving on ====================
//            REQUIRE( tLclCoord(0) =  0.2 );
//            REQUIRE( tLclCoord(1) = -0.2 );
//            REQUIRE( tNewNodeCoords(0,0) = 0.6 );
//            REQUIRE( tNewNodeCoords(1,1) = 0.6 );
//            //=======================================================================
//            Matrix< moris::IndexMat > tEdgeIndices00{{ 0,1 }};
//            std::shared_ptr< xtk::Topology > tTop00 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices00 );
//
//
//            Matrix< moris::IndexMat > tEdgeIndices01{{ 0,3 }};
//            std::shared_ptr< xtk::Topology > tTop01 = std::make_shared< xtk::Edge_Topology >( tEdgeIndices01 );
//            moris::Cell< xtk::Topology* > tTopCell(2);
//            tTopCell(0) = tTop00.get();          tTopCell(1) = tTop01.get();
//
//            moris::Cell< moris_index > tIndices(2);
//            tIndices(0) = 4;   tIndices(1) = 5;
//            moris::Cell< Matrix<DDRMat> > tLocalCoords(2);
//            tLocalCoords(0) = {{0.2}};  tLocalCoords(1) = {{-0.2}};
//            moris::Cell< Matrix<DDRMat> > tGlbCoords(2);
//            tGlbCoords(0) = tNewNodeCoords.get_row(0); tGlbCoords(1) = tNewNodeCoords.get_row(1);
//
//            tGENGeometryEngine.create_new_node_geometry_objects( tIndices, true, tTopCell, tLocalCoords, tGlbCoords );
//
//            Matrix< DDRMat > tAllCoords = {    { 0.0, 0.0 },
//                                               { 1.0, 0.0 },
//                                               { 1.0, 1.0 },
//                                               { 0.0, 1.0 },
//                                           { tRadius, 0.0 },
//                                           { 0.0, tRadius }};
//
//            Matrix< IndexMat > tInd{{ tIndices(0),tIndices(1) }};
//            tGENGeometryEngine.compute_interface_sensitivity( tInd, tAllCoords, (moris_index)0 );
//
//            GEN_Geometry_Object tGeomObj01 = tGENGeometryEngine.get_geometry_object( tIndices(0) );
//            moris::Matrix< moris::DDRMat >tSens01 = tGeomObj01.get_sensitivity_dx_dp();
//    //print( tSens01,"tSens01" );
//
//            GEN_Geometry_Object tGeomObj02 = tGENGeometryEngine.get_geometry_object( tIndices(1) );
//            moris::Matrix< moris::DDRMat >tSens02 = tGeomObj02.get_sensitivity_dx_dp();
//    //print( tSens02,"tSens02" );
//            //================================== end ===============================================
//
//            //------------------------------------------------------------------------------
    //        std::string tOutputFile = "./output_sensitivityCheck_01.exo";
    //        tInterpMesh->create_output_mesh(tOutputFile);
            delete tInterpMesh;
            delete tIntegMesh;
//            //------------------------------------------------------------------------------
        }
    }

    }   // end ge namespace
}       //end moris namespace



