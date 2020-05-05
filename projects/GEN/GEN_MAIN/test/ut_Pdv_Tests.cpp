/*
 * UT_Pdv_Checks.cpp
 *
 *  Created on: Jan 17, 2020
 *      Author: sonne
 */
#include "catch.hpp"

#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "op_equal_equal.hpp"

// GE include -----------------------------------
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Pdv_Host_Manager.hpp"
#include "cl_GEN_Dv_Enums.hpp"
#include "cl_GEN_Enums.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Phase_Table.hpp"
#include "cl_GEN_Plane.hpp"
#include "cl_GEN_Property.hpp"
#include "cl_GEN_Sphere.hpp"

// HMR includes ---------------------------------
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"

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

#include "cl_PRM_HMR_Parameters.hpp"

//------------------------------------------------------------------------------

namespace moris
{
namespace ge
{
Matrix< DDRMat > tConstValFunction( moris::Cell< Matrix< DDRMat > > & aCoeff )
            {
                return aCoeff( 0 );
            }
//------------------------------------------------------------------------------
TEST_CASE("unit test for globally consistent pdv type list","[GE],[global_pdv_type_list_check_parallel]")
{
    if(par_size() == 2)
    {
        uint tLagrangeMeshIndex = 0;
        //  HMR Parameters setup
        moris::ParameterList tParameters = prm::create_hmr_parameter_list();

        tParameters.set( "number_of_elements_per_dimension", std::string("2, 1") );
        tParameters.set( "domain_dimensions",                std::string("2, 1") );
        tParameters.set( "domain_offset",                    std::string("0, 0") );

        //        tParameters.set( "domain_sidesets", "1, 2, 3, 4" );

        tParameters.set( "truncate_bsplines", 1 );
        tParameters.set( "lagrange_orders", std::string("1") );
        tParameters.set( "lagrange_pattern", std::string("0") );
        tParameters.set( "bspline_orders", std::string("1") );
        tParameters.set( "bspline_pattern", std::string("0") );

        tParameters.set( "lagrange_output_meshes", std::string("0") );
        tParameters.set( "lagrange_input_meshes",std::string( "0") );

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

        hmr::Interpolation_Mesh_HMR *      tInterpMesh      = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );
        moris::hmr::Integration_Mesh_HMR * tIntegrationMesh = tHMR.create_integration_mesh( 1, 0, *tInterpMesh );

        mtk::Mesh_Manager tMeshManager;

        uint tHMRMeshIndex = tMeshManager.register_mesh_pair( tInterpMesh, tIntegrationMesh );
        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------
        Cell< enum GEN_DV > tPdvList0(1);
        tPdvList0(0) = GEN_DV::DENSITY0;

        std::shared_ptr< GEN_Property > tConstDensityProp0 = std::make_shared< GEN_Property >();
        tConstDensityProp0->set_parameters( { {{ 1234 }} } );
        tConstDensityProp0->set_val_function( tConstValFunction );

        Cell< std::shared_ptr< GEN_Property > > tPropertyList0(1);
        tPropertyList0(0) = tConstDensityProp0;
        //------------------------------------------------------------------------------
        Cell< enum GEN_DV > tPdvList1(1);
        tPdvList1(0) = GEN_DV::DENSITY1;

        std::shared_ptr< GEN_Property > tConstDensityProp1 = std::make_shared< GEN_Property >();
        tConstDensityProp1->set_parameters( { {{ 4321 }} } );
        tConstDensityProp1->set_val_function( tConstValFunction );

        Cell< std::shared_ptr< GEN_Property > > tPropertyList1(1);
        tPropertyList1(0) = tConstDensityProp1;
        //------------------------------------------------------------------------------

        Cell<std::shared_ptr<moris::ge::Geometry_Analytic>> tGeometry(0);
        moris::ge::Phase_Table tPhaseTable (1, moris::ge::Phase_Table_Structure::EXP_BASE_2);
        moris::ge::Geometry_Engine  tGeometryEngine(tGeometry, tPhaseTable, 2);

        //------------------------------------------------------------------------------
        tGeometryEngine.register_mesh( &tMeshManager );

        if( par_rank()==0 )
        {
            tGeometryEngine.set_pdv_types( tPdvList0 );
            tGeometryEngine.initialize_interp_pdv_host_list(  );

            tGeometryEngine.assign_ip_hosts_by_set_index( 0, tPropertyList0(0), tPdvList0(0), tHMRMeshIndex );
        }
        else if ( par_rank()==1 )
        {
            tGeometryEngine.set_pdv_types( tPdvList1 );
            tGeometryEngine.initialize_interp_pdv_host_list(  );

            tGeometryEngine.assign_ip_hosts_by_set_index( 0, tPropertyList1(0), tPdvList1(0), tHMRMeshIndex );
        }

        // ----- check the global consistent lists -----
        moris::Cell< enum GEN_DV > tCheckCell = { {GEN_DV::DENSITY0},
                                                  {GEN_DV::DENSITY1} };
        if( par_rank()==0 )
        {
            moris::Cell< enum GEN_DV > tPdvTypeList0 = tGeometryEngine.get_pdv_host_manager()->get_ip_pdv_type_list();

            REQUIRE( tPdvTypeList0(0) == tCheckCell(0) );
            REQUIRE( tPdvTypeList0(1) == tCheckCell(1) );
        }

        if( par_rank()==1 )
        {
            moris::Cell< enum GEN_DV > tPdvTypeList1 = tGeometryEngine.get_pdv_host_manager()->get_ip_pdv_type_list();

            REQUIRE( tPdvTypeList1(0) == tCheckCell(0) );
            REQUIRE( tPdvTypeList1(1) == tCheckCell(1) );
        }

    }   // end par size statement
}
//------------------------------------------------------------------------------
TEST_CASE("unit test for globally consistent pdv type list with geometry","[GE],[global_pdv_type_list_check_parallel_with_geometry]")
{
    /*
     * similar to above test but now we also add the geometry design variables on the integration nodes
     */
    if(par_size() == 2)
    {
        // TODO: implement test
    }
}


}   // end ge namespace
}       // end moris namespace
