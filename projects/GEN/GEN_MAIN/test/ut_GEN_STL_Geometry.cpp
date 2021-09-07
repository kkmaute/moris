#include "catch.hpp"
#include "cl_SDF_Generator.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "paths.hpp"

namespace moris
{
    TEST_CASE("SDF Teapot", "[SDF_Teapot]")
    {
       
        ParameterList tParameterlist = prm::create_hmr_parameter_list();
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


        hmr::HMR tHMR( tParameterlist );

        std::shared_ptr< moris::hmr::Mesh > tMesh = tHMR.create_mesh( 0 );

        for( uint k=0; k<2; ++k )
        {
            // matrices with surface element IDs
            Matrix< IndexMat > tSurfaceElements;
            tSdfGen.raycast( tMesh, tSurfaceElements );

            // get number of surface elements
            uint tNumberOfSurfaceElements = tSurfaceElements.length();

            // loop over all elements
            for( uint e=0; e<tNumberOfSurfaceElements; ++e )
            {
                // manually flag element
                tHMR.flag_element( tSurfaceElements( e ) );
            }

            tHMR.perform_refinement_based_on_working_pattern( 0  );
        }


        // calculate T-Matrices etc
        tHMR.finalize();

        // calculate SDF
        auto tField = tMesh->create_field( "SDF",  0 );

        //------------------------------------------------------------------------------

        tSdfGen.calculate_sdf( tMesh, tField->get_node_values() );

        tHMR.save_to_exodus(0, "SDF.exo",0.0 );


    }


}
