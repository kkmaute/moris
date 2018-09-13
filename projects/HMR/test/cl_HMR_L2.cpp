#include <catch.hpp>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src
#include "cl_HMR.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"


TEST_CASE("HMR_L2_Test", "[moris],[mesh],[hmr]")
{
    std::cout << "HMR_L2_Test is deactivated" << std::endl;

    /*
    // this test only works in serial so far
    if(  moris::par_size() == 1 )
        {
//-------------------------------------------------------------------------------

        SECTION("HMR Least Squares")
        {
            for( moris::uint tDimension=2; tDimension<=3; ++tDimension )
            {
                for( moris::uint tOrder=1; tOrder<=3; tOrder++ )
                {
//------------------------------------------------------------------------------
//  HMR Parameters setup
//------------------------------------------------------------------------------

                    // The parameter object controls the behavior of HMR.
                    moris::hmr::Parameters tParameters;

                    // 2D case
                    if( tDimension == 2 )
                    {
                        // We create a 2-Dimensional mesh with 2x2 elements ...
                        tParameters.set_number_of_elements_per_dimension( 2, 2 );
                    }
                    else
                    {
                        // We create a 3-Dimensional mesh with 2x2x2 elements ...
                        tParameters.set_number_of_elements_per_dimension( 2, 2, 2 );
                    }

                    // set default mesh order to tOrder
                    tParameters.set_mesh_order( tOrder );

                    // make mesh output silent
                    tParameters.set_verbose( false );

                    // B-Spline truncation is turned on by default.
                    // It is recommended to leave this setting as is.
                    tParameters.set_bspline_truncation( true );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

                    // create the HMR object by passing the settings to the constructor
                    moris::hmr::HMR tHMR( tParameters );

                    // select input pattern
                    tHMR.set_activation_pattern( tParameters.get_input_pattern() );

                    // refine the first element three times
                    // refine for three levels
                    for( uint tLevel = 0; tLevel < 3; ++tLevel )
                    {
                        tHMR.flag_element( 0 );
                        tHMR.perform_refinement();
                    }

                    // select output pattern
                    tHMR.set_activation_pattern( tParameters.get_output_pattern() );

                    // refine the last element three times
                    for( uint tLevel = 0; tLevel < 3; ++tLevel )
                    {
                        tHMR.flag_element(  tHMR.get_number_of_elements_on_proc()-1 );
                        // perform refinement
                        tHMR.perform_refinement();
                    }

//------------------------------------------------------------------------------
//  Fields
//------------------------------------------------------------------------------

                    // Create fields for all three patterns
                    auto tField0 = tHMR.create_field(
                            "Field",
                            tParameters.get_input_pattern() );
                    tField0->evaluate_function( moris::norm );

                    // perform L2 projector for output mesh
                    auto tField1 = tHMR.map_field_to_output_mesh( tField0 );

                    // calculate exact solution for reference
                    auto tExact = tHMR.create_field(
                            "Exact",
                            tField1->get_lagrange_index() );

                    tExact->evaluate_function( moris::norm );

//------------------------------------------------------------------------------
//   Test error
//------------------------------------------------------------------------------

                    // determine coefficient of determination
                    moris::real tR2 = moris::r2(
                            tExact->get_data(),
                            tField1->get_data() );

                    // perform test
                    if( tOrder == 1 )
                    {
                        REQUIRE( tR2 > 0.97 );
                    }
                    else
                    {
                        REQUIRE( tR2 > 0.99 );
                    }

//------------------------------------------------------------------------------
 //    Output
//------------------------------------------------------------------------------

                    // uncommented, since not needed for test
                    // tHMR.save_to_exodus( "Mesh.exo" );
                }
            }
        }

//-------------------------------------------------------------------------------
        }*/
}
