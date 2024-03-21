/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_STL_Geometry.cpp
 *
 */

#include "catch.hpp"
#include "cl_SDF_Generator.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "paths.hpp"

namespace moris
{
    TEST_CASE("SDF Teapot", "[SDF_Teapot]")
    {
        Parameter_List tParameterlist = prm::create_hmr_parameter_list();
        tParameterlist = prm::create_hmr_parameter_list();

        tParameterlist.set( "number_of_elements_per_dimension", "10, 10, 10" );
        tParameterlist.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
        tParameterlist.set( "domain_offset",                    "-4.9, 3.25, -1.7" );

        tParameterlist.set( "lagrange_output_meshes",           "0");

        tParameterlist.set( "lagrange_orders",  "1" );
        tParameterlist.set( "lagrange_pattern", std::string( "0" )  );
        tParameterlist.set( "bspline_orders",   "1" );
        tParameterlist.set( "bspline_pattern",  std::string( "0" )  );

        tParameterlist.set( "lagrange_to_bspline", "0" );

        tParameterlist.set( "truncate_bsplines",  1 );
        tParameterlist.set( "refinement_buffer",  1 );
        tParameterlist.set( "staircase_buffer",   1 );
        tParameterlist.set( "initial_refinement", "0");
        tParameterlist.set( "initial_refinement_pattern", "0" );

        tParameterlist.set( "use_number_aura", 1);

        tParameterlist.set( "use_multigrid",  0 );
        tParameterlist.set( "severity_level", 0 );

        std::string tObjectPath = get_base_moris_dir() + "/projects/HMR/tutorials/bracket.obj";

        // create SDF generator
        sdf::SDF_Generator tSdfGen( tObjectPath );

        // create instance of HMR
        std::shared_ptr<hmr::HMR> tHMR = std::make_shared<hmr::HMR>( tParameterlist );

        // initialize a mesh manager
        std::shared_ptr<mtk::Mesh_Manager> tMeshManager = std::make_shared<mtk::Mesh_Manager>();

        // perform HMR
        tMeshManager->set_performer(tHMR);
        tHMR->set_performer(tMeshManager);
        tHMR->perform_initial_refinement();
        tHMR->perform();

        for( uint k=0; k<2; ++k )
        {
            // matrices with surface element IDs
            Matrix< IndexMat > tSurfaceElements;
            tSdfGen.raycast( tMeshManager->get_interpolation_mesh( 0 ), tSurfaceElements );

            // get number of surface elements
            uint tNumberOfSurfaceElements = tSurfaceElements.length();

            // loop over all elements
            for( uint e=0; e<tNumberOfSurfaceElements; ++e )
            {
                // manually flag element
                // fixme: shouldn't this be: tHMR->flag_element( tSurfaceElements( e ) );
                tHMR->flag_element( e );
            }

            tHMR->perform_refinement_based_on_working_pattern( 0  );
        }

        // calculate SDF
        auto tField = tHMR->create_field( "SDF", 0, 0 );

        //------------------------------------------------------------------------------

        tSdfGen.calculate_sdf( tMeshManager->get_interpolation_mesh( 0 ), tField->get_node_values() );
    }

}

