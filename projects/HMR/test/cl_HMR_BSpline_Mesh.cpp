#include <catch.hpp>

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src

#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src


// This test creates a simple refinement pattern and makes sure that each B-Spline
// is only generated once.


TEST_CASE("HMR_Bspline_Mesh", "[moris],[mesh],[hmr]")
{
//-------------------------------------------------------------------------------

    if(  moris::par_size() == 1  ||  moris::par_size() == 2  || moris::par_size() == 4 )
    {
//-------------------------------------------------------------------------------

        SECTION("B-Spline Mesh 2D: test basis uniqueness")
        {
            // create settings object
            moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements;
            if( moris::par_size() == 1 )
            {
                tNumberOfElements.set_size( 2, 1, 3 );

            }
            else if ( moris::par_size() == 2 )
            {
                tNumberOfElements.set_size( 2, 1, 6 );

            }
            else if ( moris::par_size() == 4 )
            {
                tNumberOfElements.set_size( 2, 1, 10 );

            }

            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_verbose( false );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // set buffer size to zero
            tParameters->set_buffer_size( 0 );

            // set max order to 3
            tParameters->set_mesh_orders_simple( 3 );

            // create factory
            moris::hmr::Factory tFactory;

            // max level to refine
            uint tMaxLevel = 3;

            // loop over several orders
            for( uint tOrder=1; tOrder<=3; tOrder++ )
            {

                // set buffer size to zero
                tParameters->set_buffer_size( tOrder );

                // set aura
                //tParameters->set_max_polynomial( tOrder );

                // create background mesh object
                moris::hmr::Background_Mesh_Base* tBackgroundMesh
                = tFactory.create_background_mesh( tParameters );

                // refine a few elements in the mesh
                for( moris::uint l=0; l<tMaxLevel; ++l  )
                {
                    auto tNumberOfElements
                    =  tBackgroundMesh->get_number_of_active_elements_on_proc();

                    // refine every other element
                    for( moris::luint k=0; k<tNumberOfElements; k += 3 )
                    {
                        // get element
                        moris::hmr::Background_Element_Base* tElement
                        = tBackgroundMesh->get_element( k );

                        // flag element for refinement
                        tElement->put_on_refinement_queue();
                    }

                    // refine mesh
                    tBackgroundMesh->perform_refinement();
                }

                // create B-Spline mesh
                moris::hmr::BSpline_Mesh_Base* tBSplineMesh
                    = tFactory.create_bspline_mesh( tParameters, tBackgroundMesh, 0, tOrder );

                // test basis uniqueness
                REQUIRE ( tBSplineMesh->test_for_double_basis() );

                // tidy up
                delete tBSplineMesh;
                delete tBackgroundMesh;
            }

            delete tParameters;

        }
//-------------------------------------------------------------------------------

        SECTION("B-Spline Mesh 3D: test basis uniqueness")
        {
            // create settings object
            moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements;

            if( moris::par_size() == 1 )
            {
                tNumberOfElements.set_size( 3, 1, 3 );

            }
            else if ( moris::par_size() == 2 )
            {
                tNumberOfElements.set_size( 3, 1, 6 );

            }
            else if ( moris::par_size() == 4 )
            {
                tNumberOfElements.set_size( 3, 1, 10 );

            }

            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_verbose( false );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // set buffer size to zero
            tParameters->set_buffer_size( 0 );

            // set max order to 3
            tParameters->set_mesh_orders_simple( 3 );

            // create factory
            moris::hmr::Factory tFactory;

            // max level to refine
            uint tMaxLevel = 3;

            // loop over several orders
            for( uint tOrder=1; tOrder<=3; tOrder++ )
            {

                // set buffer size to zero
                tParameters->set_buffer_size( tOrder );

                // create background mesh object
                moris::hmr::Background_Mesh_Base* tBackgroundMesh
                    = tFactory.create_background_mesh( tParameters );

                // refine a few elements in the mesh
                for( moris::uint l=0; l<tMaxLevel; ++l  )
                {
                    auto tNumberOfElements
                    =  tBackgroundMesh->get_number_of_active_elements_on_proc();

                    // refine every other element
                    for( moris::luint k=0; k<tNumberOfElements; k += 3 )
                    {
                        // get element
                        moris::hmr::Background_Element_Base* tElement
                            = tBackgroundMesh->get_element( k );

                        // flag element for refinement
                        tElement->put_on_refinement_queue();
                    }

                    // refine mesh
                    tBackgroundMesh->perform_refinement();
                }

                // create B-Spline mesh
                moris::hmr::BSpline_Mesh_Base* tBSplineMesh
                    = tFactory.create_bspline_mesh( tParameters, tBackgroundMesh, 0, tOrder );

                // test basis uniqueness
                REQUIRE ( tBSplineMesh->test_for_double_basis() );
                //std::cout << "check " << moris::par_rank() << " " << tOrder << " " << tBSplineMesh->test_for_double_basis() << std::endl;
                // tidy up
                delete tBSplineMesh;
                delete tBackgroundMesh;
            }

            delete tParameters;

        }
    }
}
