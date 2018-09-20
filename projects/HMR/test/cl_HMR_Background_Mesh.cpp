#include <catch.hpp>

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src

#include "cl_HMR_Background_Mesh.hpp" //HMR/src
#include "cl_HMR_Background_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

TEST_CASE("HMR_Background_Mesh", "[moris],[mesh],[hmr]")
{
    SECTION( "Background mesh initialization test")
    {
        if( moris::par_size() == 1 ||  moris::par_size() == 2  || moris::par_size() == 4 )
        {
            // create settings object
            moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Mat<moris::luint> tNumberOfElements = { {4}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_verbose( false );

            // set buffer size to zero
            tParameters->set_buffer_size( 0 );

            // use simple patterns
            tParameters->set_mesh_orders_simple( 2 );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory;

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh
                = tFactory.create_background_mesh( tParameters );

            // update element table
            tBackgroundMesh->collect_active_elements();

            auto tActivePattern = tBackgroundMesh->get_activation_pattern();

            if ( moris::par_size() == 1 )
            {

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                //  tBackgroundMesh->print_level_zero();
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 45);

                // must be 10
                REQUIRE( tElement->get_domain_id() == 45 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 54 );

                // must be 38
                REQUIRE( tElement->get_domain_id() == 54 );

                // must be false
                REQUIRE( ! tElement->is_active( tActivePattern ) );

                // must be true
                REQUIRE(  tElement->is_refined( tActivePattern ) );

                // must be true
                REQUIRE(  tElement->is_padding() );


            }
            else if ( moris::par_size() == 2 )
            {
                if ( moris::par_rank() == 1)
                {
                    // print zero layer for debugging
                    // tBackgroundMesh->print_level_zero();

                    // element pointer for tests
                    moris::hmr::Background_Element_Base* tElement = nullptr;

                    // ------------------------------------------------------
                    // simple element test
                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 14 );

                    // must be
                    REQUIRE( tElement->get_domain_id() == 20 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 1
                    REQUIRE( tElement->get_owner() == 1);

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 23 );

                    // must be
                    REQUIRE( tElement->get_domain_id() == 31 );

                    // must be false
                    REQUIRE( ! tElement->is_active( tActivePattern ) );

                    // must be true
                    REQUIRE( tElement->is_refined( tActivePattern ) );

                    // must be UINT_MAX
                    // (padding elements do not belong to anybody )
                    REQUIRE( tElement->get_owner() == moris::hmr::gNoProcOwner );

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 12 );

                    // must be
                    REQUIRE( tElement->get_domain_id() == 18 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 0
                    REQUIRE( tElement->get_owner() == 0 );

                    // ------------------------------------------------------
                    // test aura from proc 1 to proc 0
                    // ------------------------------------------------------

                    // get aura elements
                    moris::Mat< moris::luint > tAuraElements
                        = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 3 );

                    REQUIRE( tAuraElements.length() == 8 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_domain_id() == 18 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_domain_id() == 19 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_domain_id() == 26 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_domain_id() == 27 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 4 ) );
                    REQUIRE( tElement->get_domain_id() == 34 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 5 ) );
                    REQUIRE( tElement->get_domain_id() == 35 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 6 ) );
                    REQUIRE( tElement->get_domain_id() == 42 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 7 ) );
                    REQUIRE( tElement->get_domain_id() == 43 );

                    // ------------------------------------------------------
                    // test aura from proc 0 to proc 1
                    // ------------------------------------------------------

                    // test inverse aura
                    moris::Mat< moris::luint > tInverseAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 3 );

                    // must be 4
                    REQUIRE( tInverseAuraElements.length() == 8 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_domain_id() == 20 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_domain_id() == 21 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_domain_id() == 28 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_domain_id() == 29 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 4 ) );
                    REQUIRE( tElement->get_domain_id() == 36 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 5 ) );
                    REQUIRE( tElement->get_domain_id() == 37 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 6 ) );
                    REQUIRE( tElement->get_domain_id() == 44 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 7 ) );
                    REQUIRE( tElement->get_domain_id() == 45 );

                }
            }
            else if ( moris::par_size() == 4 )
            {
                if ( moris::par_rank() == 2)
                {
                    // print zero layer for debugging
                    //tBackgroundMesh->print_level_zero();

                    // element pointer for tests
                    moris::hmr::Background_Element_Base* tElement = nullptr;

                    // ------------------------------------------------------
                    // simple element test
                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 14 );

                    REQUIRE( tElement->get_domain_id() == 20 );

                    REQUIRE( tElement->is_active( tActivePattern ) );

                    REQUIRE( !tElement->is_refined( tActivePattern ) );

                    REQUIRE( tElement->get_owner() == 2 );

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 12 );

                    // must be 19
                    REQUIRE( tElement->get_domain_id() == 18 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 0
                    REQUIRE( tElement->get_owner() ==  0);

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 26 );

                    REQUIRE( tElement->get_domain_id() == 36 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 3
                    REQUIRE( tElement->get_owner() == 3 );

                    // ------------------------------------------------------
                    // test aura from proc 2 to proc 1
                    // ------------------------------------------------------

                    moris::Mat< moris::luint > tAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 6 );

                    // must be 4
                    REQUIRE( tAuraElements.length() == 4);


                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_domain_id() == 34 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_domain_id() == 35 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_domain_id() == 42 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_domain_id() == 43 );

                    // ------------------------------------------------------
                    // test aura from proc 1 to proc 2
                    // ------------------------------------------------------

                    moris::Mat< moris::luint > tInverseAuraElements
                                   = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 6 );

                    // must be
                    REQUIRE( tInverseAuraElements.length() == 4 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_domain_id() == 20 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_domain_id() == 21 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_domain_id() == 28 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_domain_id() == 29 );



                }
            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
    }

    SECTION( "Background mesh refinement test")
    {
        if( moris::par_size() == 1 ||  moris::par_size() == 2  || moris::par_size() == 4 )
        {
            // create settings object
            moris::hmr::Parameters * tParameters = new moris::hmr::Parameters;

            // set number of elements
            moris::Mat<moris::luint> tNumberOfElements = { {4}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_verbose( false );

            // set buffer size to zero
            tParameters->set_buffer_size( 0 );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // use simple patterns
            tParameters->set_mesh_orders_simple( 2 );

            // create factory
            moris::hmr::Factory tFactory;

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh
                = tFactory.create_background_mesh( tParameters );

            // update element table
            tBackgroundMesh->collect_active_elements();

            if ( moris::par_size() == 1 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // since no element has been refined yet, I know that element
                // element 0 is the element with ID 18
                tElement = tBackgroundMesh->get_element( 0 );

                // refine it
                tBackgroundMesh->refine_element( tElement );

                // since the list of active elements has not been updated yet,
                // In know that 5 is the element with ID 27
                tElement = tBackgroundMesh->get_element( 5 );

                // refine it
                tBackgroundMesh->refine_element( tElement );

                // update element table
                tBackgroundMesh->collect_active_elements();

                // list of active elements on proc
                moris::Mat< moris::luint> tElementIDs;

                tBackgroundMesh->get_active_elements_on_proc( tElementIDs );

                REQUIRE( tElementIDs.length()  == 22 );

                REQUIRE( tElementIDs(  0 ) == 132 );
                REQUIRE( tElementIDs(  1 ) == 133 );
                REQUIRE( tElementIDs(  2 ) == 148 );
                REQUIRE( tElementIDs(  3 ) == 149 );
                REQUIRE( tElementIDs(  4 ) ==  19 );
                REQUIRE( tElementIDs(  5 ) ==  20 );
                REQUIRE( tElementIDs(  6 ) ==  21 );
                REQUIRE( tElementIDs(  7 ) ==  26 );
                REQUIRE( tElementIDs(  8 ) == 166 );
                REQUIRE( tElementIDs(  9 ) == 167 );
                REQUIRE( tElementIDs( 10 ) == 182 );
                REQUIRE( tElementIDs( 11 ) == 183 );
                REQUIRE( tElementIDs( 12 ) ==  28 );
                REQUIRE( tElementIDs( 13 ) ==  29 );
                REQUIRE( tElementIDs( 14 ) ==  34 );
                REQUIRE( tElementIDs( 15 ) ==  35 );
                REQUIRE( tElementIDs( 16 ) ==  36 );
                REQUIRE( tElementIDs( 17 ) ==  37 );
                REQUIRE( tElementIDs( 18 ) ==  42 );
                REQUIRE( tElementIDs( 19 ) ==  43 );
                REQUIRE( tElementIDs( 20 ) ==  44 );
                REQUIRE( tElementIDs( 21 ) ==  45 );

            }
            else if ( moris::par_size() == 2 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base * tElement = nullptr;

                // update element table
                tBackgroundMesh->collect_active_elements();

                if( moris::par_rank() == 0 )
                {
                    // pick element 19
                    tElement = tBackgroundMesh->get_element( 1 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 19 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                else
                {
                    // pick element 20
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 20 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement();
                if( moris::par_rank() == 0 )
                {
                    // pick element 150
                    tElement =  tBackgroundMesh->get_element( 3 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 150 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement();

                // list of active elements on proc
                moris::Mat< moris::luint> tElementIDs;

                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 25 );

                REQUIRE( tElementIDs(  0 ) ==  18 );
                REQUIRE( tElementIDs(  1 ) == 134 );
                REQUIRE( tElementIDs(  2 ) == 135 );
                REQUIRE( tElementIDs(  3 ) == 652 );
                REQUIRE( tElementIDs(  4 ) == 653 );
                REQUIRE( tElementIDs(  5 ) == 684 );
                REQUIRE( tElementIDs(  6 ) == 685 );
                REQUIRE( tElementIDs(  7 ) == 151 );
                REQUIRE( tElementIDs(  8 ) == 136 );
                REQUIRE( tElementIDs(  9 ) == 137 );
                REQUIRE( tElementIDs( 10 ) == 152 );
                REQUIRE( tElementIDs( 11 ) == 153 );
                REQUIRE( tElementIDs( 12 ) ==  21 );
                REQUIRE( tElementIDs( 13 ) ==  26 );
                REQUIRE( tElementIDs( 14 ) ==  27 );
                REQUIRE( tElementIDs( 15 ) ==  28 );
                REQUIRE( tElementIDs( 16 ) ==  29 );
                REQUIRE( tElementIDs( 17 ) ==  34 );
                REQUIRE( tElementIDs( 18 ) ==  35 );
                REQUIRE( tElementIDs( 19 ) ==  36 );
                REQUIRE( tElementIDs( 20 ) ==  37 );
                REQUIRE( tElementIDs( 21 ) ==  42 );
                REQUIRE( tElementIDs( 22 ) ==  43 );
                REQUIRE( tElementIDs( 23 ) ==  44 );
                REQUIRE( tElementIDs( 24 ) ==  45 );

            }
            else if ( moris::par_size() == 4 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base * tElement = nullptr;

                // update element table
                tBackgroundMesh->collect_active_elements();

                if( moris::par_rank() == 0 )
                {
                    // pick element 19
                    tElement = tBackgroundMesh->get_element( 1 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 19 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                if( moris::par_rank() == 2 )
                {
                    // pick element 20
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 20 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // call refinement procedure
                tBackgroundMesh->perform_refinement();

                if( moris::par_rank() == 0 )
                {
                    // pick element 134
                    tElement =  tBackgroundMesh->get_element( 3 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_domain_id() == 150 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // call refinement procedure
                tBackgroundMesh->perform_refinement();

                // list of active elements on proc including aura
                moris::Mat<moris::luint> tElementIDs;

                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on all four procs
                REQUIRE( tElementIDs.length() == 25 );

                REQUIRE( tElementIDs(  0 ) ==  18 );
                REQUIRE( tElementIDs(  1 ) == 134 );
                REQUIRE( tElementIDs(  2 ) == 135 );
                REQUIRE( tElementIDs(  3 ) == 652 );
                REQUIRE( tElementIDs(  4 ) == 653 );
                REQUIRE( tElementIDs(  5 ) == 684 );
                REQUIRE( tElementIDs(  6 ) == 685 );
                REQUIRE( tElementIDs(  7 ) == 151 );
                REQUIRE( tElementIDs(  8 ) == 136 );
                REQUIRE( tElementIDs(  9 ) == 137 );
                REQUIRE( tElementIDs( 10 ) == 152 );
                REQUIRE( tElementIDs( 11 ) == 153 );
                REQUIRE( tElementIDs( 12 ) ==  21 );
                REQUIRE( tElementIDs( 13 ) ==  26 );
                REQUIRE( tElementIDs( 14 ) ==  27 );
                REQUIRE( tElementIDs( 15 ) ==  28 );
                REQUIRE( tElementIDs( 16 ) ==  29 );
                REQUIRE( tElementIDs( 17 ) ==  34 );
                REQUIRE( tElementIDs( 18 ) ==  35 );
                REQUIRE( tElementIDs( 19 ) ==  36 );
                REQUIRE( tElementIDs( 20 ) ==  37 );
                REQUIRE( tElementIDs( 21 ) ==  42 );
                REQUIRE( tElementIDs( 22 ) ==  43 );
                REQUIRE( tElementIDs( 23 ) ==  44 );
                REQUIRE( tElementIDs( 24 ) ==  45 );
            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
    }
}
