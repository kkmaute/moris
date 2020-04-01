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
#include "cl_GEN_Design_Variable_Interface.hpp"
#include "cl_GEN_Dv_Enums.hpp"
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Field.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Phase_Table.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Sphere.hpp"

// HMR includes ---------------------------------
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"

// MTK includes ---------------------------------
#include "cl_Mesh_Enums.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Writer_Exodus.hpp"

// XTK include ----------------------------------
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Edge_Topology.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"

#include "cl_MSI_Design_Variable_Interface.hpp"

#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "cl_PRM_HMR_Parameters.hpp"


//------------------------------------------------------------------------------

namespace moris
{
    namespace ge
    {

    Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > >         & aCoeff )
    {
        return aCoeff( 0 );
    }

    Matrix< DDRMat > tFieldPropFunc0( moris::Cell< Matrix< DDRMat > > & aCoeff )
    {
        // fill all values with a constant: 1234
        uint tSize = aCoeff.size();
        Matrix< DDRMat > tAllVals(tSize,1);
        for(uint i=0; i<tSize; i++)
        {
            tAllVals(i) = 1234;
        }
        return tAllVals;
    }

    Matrix< DDRMat > tFieldPropFunc1( moris::Cell< Matrix< DDRMat > > & aCoeff )
    {
        // fill all values with a constant: 4321
        uint tSize = aCoeff.size();
        Matrix< DDRMat > tAllVals(tSize,1);
        for(uint i=0; i<tSize; i++)
        {
            tAllVals(i) = 4321;
        }
        return tAllVals;
    }

//------------------------------------------------------------------------------
    TEST_CASE("pdv_test_00","[GE],[pdv_check_00]")
    {
        /*
         * testing the interpolation pdv hosts in serial - parallel tests done in associated unit test
         *  - assign different density values to two materials on same mesh
         */

        if(par_size()<=1)
        {
            uint tLagrangeMeshIndex = 0;
            //  HMR Parameters setup
            moris::ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::string("2, 2, 2") );
            tParameters.set( "domain_dimensions",                std::string("2, 2, 2") );
            tParameters.set( "domain_offset",                    std::string("-1, -1, -1") );

            tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4, 5, 6") );

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

            auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >

            tHMR.perform_initial_refinement( 0 );

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();

            tHMR.finalize();

            std::shared_ptr< hmr::Interpolation_Mesh_HMR >      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );

            mtk::Mesh_Manager tMeshManager;

//            uint tHMRMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh.get(), tIntegrationMesh.get() );
            //------------------------------------------------------------------------------
            //------------------------------------------------------------------------------
            Cell< enum GEN_DV > tPdvList(3);
            tPdvList(0) = GEN_DV::DENSITY0;
            tPdvList(1) = GEN_DV::DENSITY1;
            tPdvList(2) = GEN_DV::DENSITY2;

            std::shared_ptr< GEN_Property > tConstDensityProp0 = std::make_shared< GEN_Property >();
            tConstDensityProp0->set_parameters( { {{ 1234 }} } );
            tConstDensityProp0->set_val_function( tConstValFunction );

            std::shared_ptr< GEN_Property > tConstDensityProp1 = std::make_shared< GEN_Property >();
            tConstDensityProp1->set_parameters( { {{ 4321 }} } );
            tConstDensityProp1->set_val_function( tConstValFunction );

            std::shared_ptr< GEN_Property > tConstDensityProp2 = std::make_shared< GEN_Property >();
            tConstDensityProp2->set_parameters( { {{ 1000 }} } );
            tConstDensityProp2->set_val_function( tConstValFunction );

            Cell< std::shared_ptr< GEN_Property > > tPropertyList(3);
            tPropertyList(0) = tConstDensityProp0;
            tPropertyList(1) = tConstDensityProp1;
            tPropertyList(2) = tConstDensityProp2;

            //------------------------------------------------------------------------------
            real tRadius  = 0.749;
            real tXCenter = 0.0;
            real tYCenter = 0.0;
            real tZCenter = 0.0;

            Sphere tSphere( tRadius, tXCenter, tYCenter, tZCenter );
            //------------------------------------------------------------------------------
            uint tNumDims = 3;

            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tSphere };

            moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(), Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGeometryEngine( tGeometryVector, tPhaseTable, tNumDims );

            xtk::Model tXTKModel( tNumDims, tInterpMesh.get(), &tGeometryEngine );
            tXTKModel.mVerbose = false;

            Cell<enum Subdivision_Method> tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh( );
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh( );
            //------------------------------------------------------------------------------
            bool tOutputXTKmesh = false;
            if (tOutputXTKmesh)
            {
                tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
                xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
                // Write mesh
                moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
                writer.write_mesh("", "aaaaa_pdvGeomCheck.exo");
                // Write the fields
                writer.set_time(1.0);
                writer.close_file();
            }
            //------------------------------------------------------------------------------
            uint tEnrMeshIndex = tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            tGeometryEngine.register_mesh( &tMeshManager );

            tGeometryEngine.set_pdv_types( tPdvList );
            tGeometryEngine.initialize_interp_pdv_host_list(  );

            // assign material property to the circle (density)
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_c_p0", tPropertyList(0), tPdvList(0), tEnrMeshIndex );
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_n_p0", tPropertyList(0), tPdvList(0), tEnrMeshIndex );

            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_c_p1", tPropertyList(1), tPdvList(1), tEnrMeshIndex );
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_n_p1", tPropertyList(1), tPdvList(1), tEnrMeshIndex );
            //------------------------------------------------------------------------------
            // get the vertex indices of the circle to use below
            Matrix< IndexMat > tAllVertIndices;
            moris::mtk::Set* tTempSet = tEnrIntegMesh.get_set_by_name( "HMR_dummy_c_p0" );
            moris::Cell< mtk::Cluster const * > tTempClusters = tTempSet->get_clusters_on_set();

            uint tNumClusters = tTempClusters.size();

            for(uint iClust=0; iClust<tNumClusters; iClust++)
            {
                moris::mtk::Cell const & tIPCell = tTempClusters(iClust)->get_interpolation_cell();

                moris::Cell< moris::mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();
                uint tNumVerts = tVertices.size();

                uint tOldSize = tAllVertIndices.length();

                tAllVertIndices.resize( tOldSize + tNumVerts, 1 );

                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    tAllVertIndices( tOldSize + iVert ) = tVertices( iVert )->get_index();
                }
            }
            // ---------- check the circle density values ----------
            GEN_Design_Variable_Interface tDVInterface( tGeometryEngine.get_pdv_host_manager() );

            moris::Cell< Matrix<DDRMat> > tDvVals;
            moris::Cell< moris::Matrix< DDSMat > > tIsActive;

            tDVInterface.get_ip_pdv_value( tAllVertIndices, tPdvList, tDvVals, tIsActive );

            for( uint i=0; i<tAllVertIndices.length(); i++ )
            {
                REQUIRE( tDvVals(0)(i)   == 1234 );     // tDvVal(DvType)(vertexIndex)
                REQUIRE( tIsActive(i)(0) == 1 );        // tIsActive(vertexIndex)(DvType)
            }
            // ------------------------ end ------------------------
        }
    }
//------------------------------------------------------------------------------

    TEST_CASE("pdv test 01","[GE],[pdv_check_01]")
    {
        /*
         * testing the integration node pdv functionalities in serial
         * - assign PDV types (UX, UY) to the integration nodes which are at the geometry boundaries
         * - check assignment
         */
        if(par_size()<=1)
        {
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

            // create MORIS mesh using MTK database
            mtk::MtkMeshData aMeshData( aNumElemTypes );
            aMeshData.CreateAllEdgesAndFaces = true;
            aMeshData.SpatialDim = & aNumDim;
            aMeshData.ElemConn(0) = & aElementConnQuad;
            aMeshData.NodeCoords = & aCoords;
            aMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
            aMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;

            moris::mtk::Interpolation_Mesh* tInterpMesh = moris::mtk::create_interpolation_mesh( MeshType::STK, aMeshData );
            moris::mtk::Integration_Mesh*   tIntegMesh1 = moris::mtk::create_integration_mesh(MeshType::STK,aMeshData,tInterpMesh);

            mtk::Mesh_Manager tMeshManager;
            //------------------------------------------------------------------------------
            real tRadius  = 0.5001;
            real tXCenter = 0.0;
            real tYCenter = 0.0;

            Circle tCircle( tRadius, tXCenter, tYCenter );

            uint tNumDims = 2;

            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tCircle };

            moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(), Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGeometryEngine( tGeometryVector, tPhaseTable, tNumDims );

            xtk::Model tXTKModel( tNumDims, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;

            Cell<enum Subdivision_Method> tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( EntityRank::NODE, 0 );

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh( );
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh( );

            //------------------------------------------------------------------------------
            bool tOutputXTKMesh = false;
            if (tOutputXTKMesh)
            {
                // Write mesh
                moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
                writer.write_mesh("", "aaaaa_integrationMesh.exo");
                // Write the fields
                writer.set_time(0.0);
                writer.close_file();
                /*
                 * this is how to output the mesh if not using the Exodus Writer above
                 *
                 * xtk::Output_Options tOutputOptions;
                 * tOutputOptions.mAddNodeSets = false;
                 * tOutputOptions.mAddSideSets = false;
                 * tOutputOptions.mAddClusters = false;
                 *
                 * tIntegMesh1 = tXTKModel.get_output_mesh( tOutputOptions );
                 * std::string tMeshOutputFile = "aaaaa_interpolationMesh.e";
                 *
                 * std::cout<<"-----------------------------------------"<<std::endl;
                 * std::cout<<"output mesh file name:  "<<tMeshOutputFile<<std::endl;
                 * std::cout<<"-----------------------------------------"<<std::endl;
                 *
                 * tIntegMesh1->create_output_mesh(tMeshOutputFile);
                 * delete tIntegMesh1;
                 */
            }
            //------------------------------------------------------------------------------
            uint tEnrMeshIndex = tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );

            tGeometryEngine.register_mesh( &tMeshManager );

            //------------------------------------------------------------------------------
            Cell< enum GEN_DV > tPdvList(2);
            tPdvList(0) = GEN_DV::DENSITY0;
            tPdvList(1) = GEN_DV::DENSITY1;

            std::shared_ptr< GEN_Property > tConstDensityProp0 = std::make_shared< GEN_Property >();
            tConstDensityProp0->set_parameters( { {{ 1234 }} } );
            tConstDensityProp0->set_val_function( tConstValFunction );

            std::shared_ptr< GEN_Property > tConstDensityProp1 = std::make_shared< GEN_Property >();
            tConstDensityProp1->set_parameters( { {{ 4321 }} } );
            tConstDensityProp1->set_val_function( tConstValFunction );

            Cell< std::shared_ptr< GEN_Property > > tPropertyList(2);
            tPropertyList(0) = tConstDensityProp0;
            tPropertyList(1) = tConstDensityProp1;
            //------------------------------------------------------------------------------
            tGeometryEngine.set_pdv_types( tPdvList );
            tGeometryEngine.initialize_interp_pdv_host_list(  );

            tGeometryEngine.assign_ip_hosts_by_set_index( 0, tPropertyList(0), tPdvList(0), tEnrMeshIndex );
            tGeometryEngine.assign_ip_hosts_by_set_index( 1, tPropertyList(1), tPdvList(1), tEnrMeshIndex );

            tGeometryEngine.initialize_integ_pdv_host_list(  );

            // -------- check the integration node values ----------
            GEN_Design_Variable_Interface tDVInterface( tGeometryEngine.get_pdv_host_manager() );

            moris::Cell< Matrix<DDRMat> > tDvVals;
            moris::Cell< moris::Matrix< DDSMat > > tIsActive;

            Matrix< IndexMat > tIntegVerts{{10, 11, 12}};
            Cell< enum GEN_DV > tList(2);
            tList(0) = GEN_DV::XCOORD;  tList(1) = GEN_DV::YCOORD;

            tDVInterface.get_ig_pdv_value( tIntegVerts, tList, tDvVals, tIsActive );

            REQUIRE( tDvVals(0)(0) == Approx(0.5001) );         // tDvVal(DvType)(vertexIndex)
            REQUIRE( tDvVals(1)(0) == Approx(0.0000) );
            REQUIRE( tIsActive(0)(0) == 1 );                    // tIsActive(vertexIndex)(DvType)
            REQUIRE( tIsActive(0)(1) == 1 );

            REQUIRE( tDvVals(0)(1) == Approx(0.292952) );
            REQUIRE( tDvVals(1)(1) == Approx(0.292952) );
            REQUIRE( tIsActive(1)(0) == 1 );
            REQUIRE( tIsActive(1)(1) == 1 );

            REQUIRE( tDvVals(0)(2) == Approx(0.0000) );
            REQUIRE( tDvVals(1)(2) == Approx(0.5001) );
            REQUIRE( tIsActive(2)(0) == 1 );
            REQUIRE( tIsActive(2)(1) == 1 );
            // ------------------------ end ------------------------
        }   // end par_size() <= 1
    }       // end test case
//------------------------------------------------------------------------------
    TEST_CASE("pdv_test_02","[GE],[pdv_check_02]")
    {
        /*
         * testing the interpolation pdv hosts in serial
         *  - assign different density values to two materials on same mesh, one via a GEN_Field class and the other via a property pointer
         *  - check the IP "changing" list from the interface to be sure the PDV from the field is not on it
         *  - create the IG PDVs and set XCOORD to be "unchanging"
         *  - check the IG "changing" list from the interface to be sure XCOORD is not on it
         */

        if(par_size()<=1)
        {
            uint tLagrangeMeshIndex = 0;
            //  HMR Parameters setup
            moris::ParameterList tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", std::string("2, 2, 2") );
            tParameters.set( "domain_dimensions",                std::string("2, 2, 2") );
            tParameters.set( "domain_offset",                    std::string("-1, -1, -1") );

            tParameters.set( "domain_sidesets", std::string("1, 2, 3, 4, 5, 6") );

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

            auto tDatabase = tHMR.get_database(); // std::shared_ptr< Database >

            tHMR.perform_initial_refinement( 0 );

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();

            tHMR.finalize();

            std::shared_ptr< hmr::Interpolation_Mesh_HMR >      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
            std::shared_ptr< moris::hmr::Integration_Mesh_HMR > tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );

            mtk::Mesh_Manager tMeshManager;
//            uint tHMRMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh.get(), tIntegrationMesh.get() );
            //------------------------------------------------------------------------------
            Cell< enum GEN_DV > tPdvList(2);
            tPdvList(0) = GEN_DV::DENSITY0;
            tPdvList(1) = GEN_DV::DENSITY1;

            std::shared_ptr< GEN_Property > tConstDensityProp0 = std::make_shared< GEN_Property >();
            tConstDensityProp0->set_val_function( tFieldPropFunc0 );

            std::shared_ptr< GEN_Property > tConstDensityProp1 = std::make_shared< GEN_Property >();
            tConstDensityProp1->set_parameters( { {{ 4321 }} } );
            tConstDensityProp1->set_val_function( tConstValFunction );

            Cell< std::shared_ptr< GEN_Property > > tPropertyList(2);
            tPropertyList(0) = tConstDensityProp0;
            tPropertyList(1) = tConstDensityProp1;
            //------------------------------------------------------------------------------
            real tRadius  = 0.749;
            real tXCenter = 0.0;
            real tYCenter = 0.0;
            real tZCenter = 0.0;

            Sphere tSphere( tRadius, tXCenter, tYCenter, tZCenter );
            //------------------------------------------------------------------------------
            uint tNumDims = 3;

            moris::Cell<moris::ge::GEN_Geometry*> tGeometryVector = { &tSphere };

            moris::ge::GEN_Phase_Table      tPhaseTable( tGeometryVector.size(), Phase_Table_Structure::EXP_BASE_2 );
            moris::ge::GEN_Geometry_Engine  tGeometryEngine( tGeometryVector, tPhaseTable, tNumDims );

            xtk::Model tXTKModel( tNumDims, tInterpMesh.get(), &tGeometryEngine );
            tXTKModel.mVerbose = false;

            Cell<enum Subdivision_Method> tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            tXTKModel.decompose(tDecompositionMethods);

            tXTKModel.perform_basis_enrichment( EntityRank::BSPLINE, 0 );

            xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh( );
            xtk::Enriched_Integration_Mesh   & tEnrIntegMesh  = tXTKModel.get_enriched_integ_mesh( );
            //------------------------------------------------------------------------------
            bool tOutputXTKmesh = true;
            if (tOutputXTKmesh)
            {
                tXTKModel.perform_basis_enrichment(EntityRank::BSPLINE,0);
                xtk::Enriched_Interpolation_Mesh & tEnrInterpMesh = tXTKModel.get_enriched_interp_mesh();
                xtk::Enriched_Integration_Mesh   & tEnrIntegMesh = tXTKModel.get_enriched_integ_mesh();
                // Write mesh
                moris::mtk::Writer_Exodus writer(&tEnrIntegMesh);
                writer.write_mesh("", "aaaaa_pdvGeomCheck.exo");
                // Write the fields
                writer.set_time(1.0);
                writer.close_file();
            }
            //------------------------------------------------------------------------------
            uint tEnrMeshIndex = tMeshManager.register_mesh_pair( &tEnrInterpMesh, &tEnrIntegMesh );
            //------------------------------------------------------------------------------
            // create field for DENSITY0
            std::shared_ptr< GEN_Field > tField0 = std::make_shared< GEN_Field >( tPropertyList(0), &tMeshManager );
            //------------------------------------------------------------------------------
            tGeometryEngine.register_mesh( &tMeshManager );

            tGeometryEngine.set_pdv_types( tPdvList );
            tGeometryEngine.initialize_interp_pdv_host_list(  );

            // assign material property to the circle (density)
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_c_p0", tField0, tPdvList(0), tEnrMeshIndex );
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_n_p0", tField0, tPdvList(0), tEnrMeshIndex );

            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_c_p1", tPropertyList(1), tPdvList(1), tEnrMeshIndex );
            tGeometryEngine.assign_ip_hosts_by_set_name( "HMR_dummy_n_p1", tPropertyList(1), tPdvList(1), tEnrMeshIndex );
            //------------------------------------------------------------------------------
            // create the IG hosts, set XCOORD to unchanging
            tGeometryEngine.initialize_integ_pdv_host_list(  );

            tGeometryEngine.mark_ig_dv_type_as_unchanging( GEN_DV::XCOORD );
            //------------------------------------------------------------------------------
            // get the vertex indices of the circle to use below
            Matrix< IndexMat > tAllVertIndices;
            moris::mtk::Set* tTempSet = tEnrIntegMesh.get_set_by_name( "HMR_dummy_c_p0" );
            moris::Cell< mtk::Cluster const * > tTempClusters = tTempSet->get_clusters_on_set();

            uint tNumClusters = tTempClusters.size();

            for(uint iClust=0; iClust<tNumClusters; iClust++)
            {
                moris::mtk::Cell const & tIPCell = tTempClusters(iClust)->get_interpolation_cell();

                moris::Cell< moris::mtk::Vertex * > tVertices = tIPCell.get_vertex_pointers();
                uint tNumVerts = tVertices.size();

                uint tOldSize = tAllVertIndices.length();

                tAllVertIndices.resize( tOldSize + tNumVerts, 1 );

                for(uint iVert=0; iVert<tNumVerts; iVert++)
                {
                    tAllVertIndices( tOldSize + iVert ) = tVertices( iVert )->get_index();
                }
            }
            // ---------- check the circle density values ----------
            GEN_Design_Variable_Interface tDVInterface( tGeometryEngine.get_pdv_host_manager() );

            moris::Cell< Matrix<DDRMat> > tDvVals;
            moris::Cell< moris::Matrix< DDSMat > > tIsActive;

            tDVInterface.get_ip_pdv_value( tAllVertIndices, tPdvList, tDvVals, tIsActive );

            for( uint i=0; i<tAllVertIndices.length(); i++ )
            {
                REQUIRE( tDvVals(0)(i)   == 1234 );     // tDvVal(DvType)(vertexIndex)
                REQUIRE( tIsActive(i)(0) == 1 );        // tIsActive(vertexIndex)(DvType)
            }
            // ------------------------ end ------------------------

            // ------- check the IP list of changing types ---------
            Cell< enum GEN_DV > tIPChangingList;
            tDVInterface.get_ip_requested_dv_types( tIPChangingList );

            // check that the changing type is DENSITY1
            for( uint i=0; i<tIPChangingList.size(); i++ )
            {
                REQUIRE( tIPChangingList(i) != GEN_DV::DENSITY0 );
            }
            // ------------------------ end ------------------------

            // -------- check the IG list of changing types --------
            Cell< enum GEN_DV > tIGChangingList;
            tDVInterface.get_ig_requested_dv_types( tIGChangingList );

            // check that the changing types are YCOORD and ZCOORD
            for( uint i=0; i<tIGChangingList.size(); i++ )
            {
                REQUIRE( tIGChangingList(i) != GEN_DV::XCOORD );
            }
            // ------------------------ end ------------------------
        }
    }
//------------------------------------------------------------------------------
    }   // end ge namespace
}       //end moris namespace



