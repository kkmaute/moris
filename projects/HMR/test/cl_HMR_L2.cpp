#include <catch.hpp>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Matrix.hpp" // LNA/src
#include "cl_HMR.hpp"
#include "fn_r2.hpp"
#include "fn_norm.hpp"

moris::real
LevelSetFunction( const moris::Mat< moris::real > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

TEST_CASE("HMR_L2_Test", "[moris],[mesh],[hmr]")
{
    // can only perform test for 1, 2 or 4 procs
    if( moris::par_size() == 1 || moris::par_size() == 2 || moris::par_size() == 4 )
    {
        // do this test for 2 and 3 dimensions
        for( moris::uint tDimension=2; tDimension<=3; ++tDimension )
        {

            // do this for first, second and third order
            for( moris::uint tOrder=1; tOrder<=3; tOrder++ )
            {
//------------------------------------------------------------------------------
//  HMR Parameters setup
//------------------------------------------------------------------------------

                // The parameter object controls the behavior of HMR.
                moris::hmr::Parameters tParameters;

                moris::Mat< moris::luint > tNumberOfElements;

                // set element size
                if( moris::par_size() == 1 )
                {
                    tNumberOfElements.set_size( tDimension, 1, 2 );
                }
                else if ( moris::par_size() == 2 )
                {
                    tNumberOfElements.set_size( tDimension, 1, 6 );
                }
                else if ( moris::par_size() == 4 )
                {
                    tNumberOfElements.set_size( tDimension, 1, 10 );
                }

                // set values to parameters
                tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

                // make mesh output silent
                tParameters.set_verbose( false );

                // B-Spline truncation is turned on by default.
                // It is recommended to leave this setting as is.
                tParameters.set_bspline_truncation( true );

                // set mesh order
                tParameters.set_mesh_order( tOrder );

//------------------------------------------------------------------------------
//  HMR Initialization
//------------------------------------------------------------------------------

                // create the HMR object by passing the settings to the constructor
                moris::hmr::HMR tHMR( tParameters );

                // select input pattern
                tHMR.set_activation_pattern( tParameters.get_input_pattern() );


                // manually select output pattern
                tHMR.mBackgroundMesh->set_activation_pattern( tParameters.get_input_pattern() );


                // refine the first element three times
                // fixme: change this to 3
                for( uint tLevel = 0; tLevel < 1; ++tLevel )
                {
                    tHMR.flag_element( 0 );
                    tHMR.mBackgroundMesh->perform_refinement();
                }

                // manually select output pattern
                tHMR.mBackgroundMesh->set_activation_pattern( tParameters.get_output_pattern() );

                // refine the last element three times
                // fixme: change this to 3
                for( uint tLevel = 0; tLevel < 1; ++tLevel )
                {
                    tHMR.flag_element(  tHMR.get_number_of_elements_on_proc()-1 );
                    // perform refinement
                    tHMR.mBackgroundMesh->perform_refinement();
                }

                // manually create union
                tHMR.unite_patterns(
                        tParameters.get_input_pattern(),
                        tParameters.get_output_pattern(),
                        tParameters.get_union_pattern() );

                // update background mesh
                // test if max polynomial is 3
                if ( tParameters.get_max_polynomial() > 2 )
                {
                    // activate extra pattern for exodus
                    tHMR.add_extra_refinement_step_for_exodus();
                }

                //tHMR.mBackgroundMesh->save_to_vtk("Background.vtk");
                //tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

                tHMR.update_meshes();
//------------------------------------------------------------------------------
//  Fields
//------------------------------------------------------------------------------

                // create pointer to input mesh
                moris::hmr::Mesh * tInputMesh = tHMR.create_input_mesh();

                // create pointer to input field
                moris::mtk::Field * tInputField
                    = tInputMesh->create_field( "LevelSet" );

                // evaluate function
                tInputField->evaluate_scalar_function( LevelSetFunction );

                // create pointer to output mesh
                moris::hmr::Mesh * tOutputMesh = tHMR.create_output_mesh();

                // calculate exact value
                moris::mtk::Field * tExact = tOutputMesh->create_field( "Exact" );

                tExact->evaluate_scalar_function( LevelSetFunction );

                // map input to output
                moris::mtk::Field * tOutputField
                    = tHMR.map_field_on_mesh( tInputField, tOutputMesh );

//------------------------------------------------------------------------------
//   Test error
//------------------------------------------------------------------------------

                // determine coefficient of determination
                moris::real tR2 = moris::r2(
                        tExact->get_node_values(),
                        tOutputField->get_node_values() );

                // std::cout << "R2 " << tR2 << std::endl;

                // perform test
                if( tOrder == 1 )
                {
                     REQUIRE( tR2 > 0.97 );
                }
                else
                {
                    REQUIRE( tR2 > 0.99 );
                }

                // delete input field pointer
                delete tInputField;

                // delete output field
                delete tOutputField;

                // delete exact field
                delete tExact;

                // delete input mesh pointer
                delete tInputMesh;

                // delete output mesh pointer
                delete tOutputMesh;

            } // end order loop
        } // end dimension loop
    } // end parallel
}
