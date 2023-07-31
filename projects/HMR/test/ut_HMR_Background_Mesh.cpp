/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_Background_Mesh.cpp
 *
 */

#include "cl_HMR_Background_Mesh.hpp" //HMR/src

#include <catch.hpp>
#include "cl_HMR_Background_Element.hpp" //HMR/src
#include "cl_HMR_Factory.hpp" //HMR/src
#include "cl_HMR_Parameters.hpp" //HMR/src

#include "cl_Communication_Manager.hpp" //COM/src
#include "cl_Communication_Tools.hpp" //COM/src
#include "typedefs.hpp" //COR/src
#include "cl_Matrix.hpp" //LINALG/src
#include "cl_HMR_Background_Mesh.hpp" //HMR/src

#include "cl_HMR.hpp" //HMR/src
#include "cl_HMR_Database.hpp" //HMR/src

TEST_CASE("HMR_Background_Mesh", "[moris],[mesh],[hmr],[Background_Mesh],[Background_Mesh1]")
{
    SECTION( "Background mesh initialization test")
                            {
        if( moris::par_size() == 1 ||  moris::par_size() == 2  || moris::par_size() == 4 )
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 1;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_severity_level( 0 );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 1 );
            tParameters->set_staircase_buffer( 1 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            //tBackgroundMesh->collect_active_elements();

            uint tActivePattern = tBackgroundMesh->get_activation_pattern();

            //            tBackgroundMesh->save_to_vtk( "BackgorundMesh_2Proc.vtk");

            if ( moris::par_size() == 1 )
            {
                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                //  tBackgroundMesh->print_level_zero();
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 45 );

                // must be 10
                REQUIRE( tElement->get_hmr_id() == 45 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 54 );

                // must be 38
                REQUIRE( tElement->get_hmr_id() == 54 );

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
                    REQUIRE( tElement->get_hmr_id() == 20 );

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
                    REQUIRE( tElement->get_hmr_id() == 31 );

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
                    REQUIRE( tElement->get_hmr_id() == 18 );

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
                    moris::Matrix< moris::DDLUMat > tAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 3 );

                    REQUIRE( tAuraElements.length() == 8 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 18 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 19 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 26 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 27 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 4 ) );
                    REQUIRE( tElement->get_hmr_id() == 34 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 5 ) );
                    REQUIRE( tElement->get_hmr_id() == 35 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 6 ) );
                    REQUIRE( tElement->get_hmr_id() == 42 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 7 ) );
                    REQUIRE( tElement->get_hmr_id() == 43 );

                    // ------------------------------------------------------
                    // test aura from proc 0 to proc 1
                    // ------------------------------------------------------

                    // test inverse aura
                    moris::Matrix< moris::DDLUMat > tInverseAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 3 );

                    // must be 4
                    REQUIRE( tInverseAuraElements.length() == 8 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 20 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 21 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 28 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 29 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 4 ) );
                    REQUIRE( tElement->get_hmr_id() == 36 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 5 ) );
                    REQUIRE( tElement->get_hmr_id() == 37 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 6 ) );
                    REQUIRE( tElement->get_hmr_id() == 44 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 7 ) );
                    REQUIRE( tElement->get_hmr_id() == 45 );
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

                    REQUIRE( tElement->get_hmr_id() == 20 );

                    REQUIRE( tElement->is_active( tActivePattern ) );

                    REQUIRE( !tElement->is_refined( tActivePattern ) );

                    REQUIRE( tElement->get_owner() == 2 );

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 12 );

                    // must be 19
                    REQUIRE( tElement->get_hmr_id() == 18 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 0
                    REQUIRE( tElement->get_owner() ==  0);

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 26 );

                    REQUIRE( tElement->get_hmr_id() == 36 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 3
                    REQUIRE( tElement->get_owner() == 3 );

                    // ------------------------------------------------------
                    // test aura from proc 2 to proc 1
                    // ------------------------------------------------------

                    moris::Matrix< moris::DDLUMat > tAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 6 );

                    // must be 4
                    REQUIRE( tAuraElements.length() == 4);

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 34 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 35 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 42 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 43 );

                    // ------------------------------------------------------
                    // test aura from proc 1 to proc 2
                    // ------------------------------------------------------

                    moris::Matrix< moris::DDLUMat > tInverseAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 6 );

                    // must be
                    REQUIRE( tInverseAuraElements.length() == 4 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 20 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 21 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 28 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 29 );
                }
            }

            tParameters->set_severity_level( 2 );

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
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 1;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {4} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            if ( moris::par_size() == 1 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // since no element has been refined yet, I know that element
                // element 0 is the element with ID 18
                tElement = tBackgroundMesh->get_element( 0 );

                // refine it
                tBackgroundMesh->refine_element( tElement, false );

                // since the list of active elements has not been updated yet,
                // In know that 5 is the element with ID 27
                tElement = tBackgroundMesh->get_element( 5 );

                // refine it
                tBackgroundMesh->refine_element( tElement, false );

                // update element table
                tBackgroundMesh->collect_active_elements();

                //                tBackgroundMesh->save_to_vtk( "BackgorundMesh_refine.vtk");

                // list of active elements on proc
                moris::Matrix< moris::DDLUMat > tElementIDs;

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
                    REQUIRE( tElement->get_hmr_id() == 19 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                else
                {
                    // pick element 20
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 20 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement( 0 );

                if( moris::par_rank() == 0 )
                {
                    // pick element 150
                    tElement =  tBackgroundMesh->get_element( 3 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 150 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement( 0 );

                // list of active elements on proc
                moris::Matrix< moris::DDLUMat > tElementIDs;

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
                    REQUIRE( tElement->get_hmr_id() == 19 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                if( moris::par_rank() == 2 )
                {
                    // pick element 20
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 20 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // call refinement procedure
                tBackgroundMesh->perform_refinement( 0 );

                if( moris::par_rank() == 0 )
                {
                    // pick element 134
                    tElement =  tBackgroundMesh->get_element( 3 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 150 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // call refinement procedure
                tBackgroundMesh->perform_refinement( 0 );

                // list of active elements on proc including aura
                moris::Matrix< moris::DDLUMat > tElementIDs;

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

TEST_CASE("HMR_Background_Mesh_Activation_Pattern", "[moris],[mesh],[hmr],[Background_Mesh_Activation_Pattern],[Background_Mesh]")
{
    if( moris::par_size() == 1 )
    {
        // create settings object
        auto tParameters = new moris::hmr::Parameters;

        // set number of elements
        tParameters->set_number_of_elements_per_dimension( { {4}, {4} } );

        // Processor Decomposition Method
        uint tDecompMethod = 1;
        tParameters->set_processor_decomp_method(tDecompMethod);

        // set buffer size to zero
        tParameters->set_refinement_buffer( 1 );
        tParameters->set_staircase_buffer( 1 );

        tParameters->set_lagrange_orders  ( { {2} });
        tParameters->set_lagrange_patterns({ {0} });

        // deactivate truncation
        tParameters->set_bspline_truncation( false );

        // create factory
        moris::hmr::Factory tFactory( tParameters );

        // create background mesh object
        moris::hmr::Background_Mesh_Base * tBackgroundMesh = tFactory.create_background_mesh();

        //----------------------------------------------------------------------------------------------------------
        // Work on activation pattern 0 mesh

        tBackgroundMesh->set_activation_pattern( 0 );

        // element 0 is the element with ID 18
        tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 0 );

        tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 0 );

        //        tBackgroundMesh->save_to_vtk( "BackgorundMesh_Procref1.vtk");

        //----------------------------------------------------------------------------------------------------------
        // Work on activation pattern 1 mesh

        tBackgroundMesh->set_activation_pattern( 1 );

        // element 0 is the element with ID 18
        tBackgroundMesh->get_element( 15 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 1 );

        tBackgroundMesh->get_element( 18 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 1 );

        //        tBackgroundMesh->save_to_vtk( "BackgorundMesh_Procref2.vtk");

        //----------------------------------------------------------------------------------------------------------
        // check elements on pattern 0

        tBackgroundMesh->set_activation_pattern( 0 );

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs_1;

        // get hmr ids of active elements
        tBackgroundMesh->get_active_elements_on_proc( tElementIDs_1 );

        REQUIRE( tElementIDs_1.length()  == 28 );

        REQUIRE( tElementIDs_1(  0 ) == 584 );        REQUIRE( tElementIDs_1(  2 ) == 616 );
        REQUIRE( tElementIDs_1(  4 ) == 133 );        REQUIRE( tElementIDs_1( 10 ) == 151 );
        REQUIRE( tElementIDs_1( 11 ) ==  20 );        REQUIRE( tElementIDs_1( 16 ) == 181 );
        REQUIRE( tElementIDs_1( 17 ) ==  27 );        REQUIRE( tElementIDs_1( 27 ) ==  45 );

        //----------------------------------------------------------------------------------------------------------
        // check elements on pattern 1

        tBackgroundMesh->set_activation_pattern( 1 );

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs_2;

        // get hmr ids of active elements
        tBackgroundMesh->get_active_elements_on_proc( tElementIDs_2 );

        REQUIRE( tElementIDs_2.length()  == 28 );

        REQUIRE( tElementIDs_2(  0 ) ==  18 );        REQUIRE( tElementIDs_2(  2 ) ==  20 );
        REQUIRE( tElementIDs_2(  4 ) ==  26 );        REQUIRE( tElementIDs_2( 10 ) ==  36 );
        REQUIRE( tElementIDs_2( 11 ) == 202 );        REQUIRE( tElementIDs_2( 16 ) ==  43 );
        REQUIRE( tElementIDs_2( 17 ) == 232 );        REQUIRE( tElementIDs_2( 27 ) == 1079 );

        // delete background mesh
        delete tBackgroundMesh;

        // delete settings object
        delete tParameters;
    }
}

TEST_CASE("HMR_Background_Mesh_Unite_Pattern", "[moris],[mesh],[hmr],[Background_Mesh_Unite_Pattern],[Background_Mesh]")
{
    if( moris::par_size() == 1 )
    {
        // create settings object
        auto tParameters = new moris::hmr::Parameters;

        // Processor Decomposition Method
        uint tDecompMethod = 1;
        tParameters->set_processor_decomp_method(tDecompMethod);

        // set number of elements
        tParameters->set_number_of_elements_per_dimension( { {4}, {4} } );

        // set buffer size to zero
        tParameters->set_refinement_buffer( 1 );
        tParameters->set_staircase_buffer( 1 );

        tParameters->set_lagrange_orders  ( { {2} });
        tParameters->set_lagrange_patterns({ {0} });

        // deactivate truncation
        tParameters->set_bspline_truncation( false );

        // create factory
        moris::hmr::Factory tFactory( tParameters );

        // create background mesh object
        moris::hmr::Background_Mesh_Base * tBackgroundMesh = tFactory.create_background_mesh();

        //----------------------------------------------------------------------------------------------------------
        // Work on activation pattern 0 mesh
        tBackgroundMesh->set_activation_pattern( 0 );

        // element 0 is the element with ID 18
        tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 0 );

        tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 0 );

        //        tBackgroundMesh->save_to_vtk( "BackgorundMesh_Procref1.vtk");

        //----------------------------------------------------------------------------------------------------------
        // Work on activation pattern 1 mesh
        tBackgroundMesh->set_activation_pattern( 1 );

        // element 0 is the element with ID 18
        tBackgroundMesh->get_element( 15 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 1 );

        tBackgroundMesh->get_element( 18 )->put_on_refinement_queue();
        tBackgroundMesh->perform_refinement( 1 );

        //        tBackgroundMesh->save_to_vtk( "BackgorundMesh_Procref2.vtk");

        //----------------------------------------------------------------------------------------------------------
        // unite pattern 0 and 1 on pattern 2
        tBackgroundMesh->unite_patterns( 0, 1, 2 );

        //----------------------------------------------------------------------------------------------------------
        // check pattern 2
        tBackgroundMesh->set_activation_pattern( 2 );

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs;

        // get hmr ids of active elements
        tBackgroundMesh->get_active_elements_on_proc( tElementIDs );

        REQUIRE( tElementIDs.length()  == 40 );

        REQUIRE( tElementIDs(  0 ) == 584 );        REQUIRE( tElementIDs(  2 ) == 616 );
        REQUIRE( tElementIDs(  4 ) == 133 );        REQUIRE( tElementIDs( 10 ) == 151 );
        REQUIRE( tElementIDs( 11 ) ==  20 );        REQUIRE( tElementIDs( 16 ) == 181 );
        REQUIRE( tElementIDs( 17 ) ==  27 );        REQUIRE( tElementIDs( 26 ) == 219 );
        REQUIRE( tElementIDs( 27 ) ==  42 );        REQUIRE( tElementIDs( 28 ) ==  43 );
        REQUIRE( tElementIDs( 29 ) == 232 );        REQUIRE( tElementIDs( 35 ) == 250 );
        REQUIRE( tElementIDs( 36 ) == 1046 );       REQUIRE( tElementIDs( 39 ) == 1079 );

        //----------------------------------------------------------------------------------------------------------
        // unite pattern 0 and 1 on pattern 3
        moris::Cell< uint > tSourcePattern( 2, 0 );
        tSourcePattern( 1 ) = 1;
        tBackgroundMesh->unite_patterns( tSourcePattern, 3 );

        //----------------------------------------------------------------------------------------------------------
        // check pattern 2
        tBackgroundMesh->set_activation_pattern( 3 );

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs_1;

        // get hmr ids of active elements
        tBackgroundMesh->get_active_elements_on_proc( tElementIDs_1 );

        REQUIRE( tElementIDs_1.length()  == 40 );

        REQUIRE( tElementIDs_1(  0 ) == 584 );        REQUIRE( tElementIDs_1(  2 ) == 616 );
        REQUIRE( tElementIDs_1(  4 ) == 133 );        REQUIRE( tElementIDs_1( 10 ) == 151 );
        REQUIRE( tElementIDs_1( 11 ) ==  20 );        REQUIRE( tElementIDs_1( 16 ) == 181 );
        REQUIRE( tElementIDs_1( 17 ) ==  27 );        REQUIRE( tElementIDs_1( 26 ) == 219 );
        REQUIRE( tElementIDs_1( 27 ) ==  42 );        REQUIRE( tElementIDs_1( 28 ) ==  43 );
        REQUIRE( tElementIDs_1( 29 ) == 232 );        REQUIRE( tElementIDs_1( 35 ) == 250 );
        REQUIRE( tElementIDs_1( 36 ) == 1046 );       REQUIRE( tElementIDs_1( 39 ) == 1079 );

        // delete background mesh
        delete tBackgroundMesh;

        // delete settings object
        delete tParameters;
    }
}

TEST_CASE("HMR_Background_Mesh_refine", "[moris],[mesh],[hmr],[Background_Mesh_refine],[Background_Mesh]")
{
    if( moris::par_size() == 1 )
    {
        // create parameter object
        moris::hmr::Parameters tParameters;

        tParameters.set_processor_decomp_method(1);

        tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
        //        tParameters.set_verbose( true );
        tParameters.set_severity_level( 0 );
        tParameters.set_multigrid( false );
        tParameters.set_bspline_truncation( true );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        moris::Cell< moris::Matrix< moris::DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create HMR object
        moris::hmr::HMR tHMR( tParameters );

        for( moris::uint Ii = 0; Ii<2; Ii++)
        {
            moris::luint tNumActiveElements = tHMR.get_database()->get_background_mesh()->get_number_of_active_elements_on_proc();

            for( moris::luint Ik = 0; Ik<tNumActiveElements; Ik++)
            {
                // flag first element for refinement
                tHMR.flag_element( Ik );
            }
            tHMR.perform_refinement( 0 );
            tHMR.update_refinement_pattern( 0 );
        }

        tHMR.finalize();

        //tHMR.save_background_mesh_to_vtk( "BackgorundMesh_refine.vtk");

        tParameters.set_severity_level( 0 );
    }
}

TEST_CASE("HMR_Background_Mesh_refinement_buffer", "[moris],[mesh],[hmr],[Background_Mesh_refinement_buffer],[Background_Mesh]")
{
    if( moris::par_size() == 2 )
    {
        // create parameter object
        moris::hmr::Parameters tParameters;

        tParameters.set_processor_decomp_method(1);

        tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 } } );
        //        tParameters.set_verbose( true );
        tParameters.set_severity_level( 0 );
        tParameters.set_multigrid( false );
        tParameters.set_bspline_truncation( true );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_number_aura( false );

        moris::Cell< moris::Matrix< moris::DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create HMR object
        moris::hmr::HMR tHMR( tParameters );

        for( moris::uint Ii = 0; Ii<1; Ii++)
        {
            if( moris::par_rank() == 0)
            {
                tHMR.get_database()->get_background_mesh()->get_element( 1 )->put_on_refinement_queue();
            }
            if( moris::par_rank() == 1)
            {
                tHMR.get_database()->get_background_mesh()->get_element( 6 )->put_on_refinement_queue();
            }

            tHMR.perform_refinement( 0 );
            tHMR.update_refinement_pattern( 0 );
        }

        tHMR.finalize();

        //tHMR.save_background_mesh_to_vtk( "BackgorundMesh_refine.vtk");

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs;

        // get hmr ids of active elements
        tHMR.get_database()->get_background_mesh()->get_active_elements_on_proc( tElementIDs );

        REQUIRE( tElementIDs.numel()  == 26 );

        if( moris::par_rank()==0)
        {
            REQUIRE( tElementIDs(  0 )  == 62 );			REQUIRE( tElementIDs(  1 )  == 63 );
            REQUIRE( tElementIDs(  14 ) == 100 );			REQUIRE( tElementIDs(  15 ) == 101 );
            REQUIRE( tElementIDs(  16 ) == 19 );			REQUIRE( tElementIDs(  21 ) == 25 );
            REQUIRE( tElementIDs(  24 ) == 148 );			REQUIRE( tElementIDs(  25 ) == 149 );
        }

        if( moris::par_rank()==1)
        {
            REQUIRE( tElementIDs(  0 )   == 66 );			REQUIRE( tElementIDs(  1 )   == 67 );
            REQUIRE( tElementIDs(  7 )   == 102 );			REQUIRE( tElementIDs(  8 )   == 103 );
            REQUIRE( tElementIDs(  9 )   == 16 );			REQUIRE( tElementIDs(  17 )  == 129 );
            REQUIRE( tElementIDs(  23 )  == 141 );			REQUIRE( tElementIDs(  25 )  == 153 );
        }

        tParameters.set_severity_level( 0 );
    }
}

TEST_CASE("HMR_Background_Mesh_refinement_buffer_2", "[moris],[mesh],[hmr],[Background_Mesh_refinement_buffer_2],[Background_Mesh]")
{
    if( moris::par_size() == 4 )
    {
        // create parameter object
        moris::hmr::Parameters tParameters;

        tParameters.set_processor_decomp_method(1);

        tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 } } );
        //        tParameters.set_verbose( true );
        tParameters.set_severity_level( 0 );
        tParameters.set_multigrid( false );
        tParameters.set_bspline_truncation( true );

        tParameters.set_refinement_buffer( 2 );
        tParameters.set_staircase_buffer( 1 );

        tParameters.set_lagrange_orders  ( { {1} });
        tParameters.set_lagrange_patterns({ {0} });

        tParameters.set_bspline_orders   ( { {1} } );
        tParameters.set_bspline_patterns ( { {0} } );

        tParameters.set_number_aura( false );

        moris::Cell< moris::Matrix< moris::DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { {0} };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        // create HMR object
        moris::hmr::HMR tHMR( tParameters );

        for( moris::uint Ii = 0; Ii<1; Ii++)
        {
            if( moris::par_rank() == 0)
            {
                tHMR.get_database()->get_background_mesh()->get_element( 1 )->put_on_refinement_queue();
            }
            if( moris::par_rank() == 3)
            {
                tHMR.get_database()->get_background_mesh()->get_element( 3 )->put_on_refinement_queue();
            }

            tHMR.perform_refinement( 0 );
            tHMR.update_refinement_pattern( 0 );
        }

        tHMR.finalize();

        //tHMR.save_background_mesh_to_vtk( "BackgorundMesh_refine_4.vtk");

        // list of active elements on proc
        moris::Matrix< moris::DDLUMat > tElementIDs;

        // get hmr ids of active elements
        tHMR.get_database()->get_background_mesh()->get_active_elements_on_proc( tElementIDs );

        if( moris::par_rank() == 1)
        {
            REQUIRE( tElementIDs.numel()  == 13 );
        }
        else
        {
            REQUIRE( tElementIDs.numel()  == 16 );
        }

        tParameters.set_severity_level( 0 );
    }
}

TEST_CASE("HMR_Background_Mesh_Proc_Decomp_Min_Mesh_Interface", "[moris],[mesh],[hmr],[Proc_Decomp_Min_Mesh],[Background_Mesh]")
{
    SECTION( "Background mesh initialization test")
                    {
        if( moris::par_size() == 1 ||  moris::par_size() == 2  || moris::par_size() == 4 )
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 2;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {2}, {8} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_severity_level( 0 );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 1 );
            tParameters->set_staircase_buffer( 1 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            //tBackgroundMesh->collect_active_elements();

            uint tActivePattern = tBackgroundMesh->get_activation_pattern();

            //            tBackgroundMesh->save_to_vtk( "BackgorundMesh_2Proc.vtk");

            if ( moris::par_size() == 1 )
            {
                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                //  tBackgroundMesh->print_level_zero();
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 57 );

                // must be 45
                REQUIRE( tElement->get_hmr_id() == 57 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 58 );

                // must be 46
                REQUIRE( tElement->get_hmr_id() == 58 );

                // must be false
                REQUIRE( ! tElement->is_active( tActivePattern ) );

                // must be true
                REQUIRE(  tElement->is_refined( tActivePattern ) );

                // must be true
                REQUIRE(  tElement->is_padding() );

                //  tBackgroundMesh->print_level_zero();
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 14 );

                // must be 26
                REQUIRE( tElement->get_hmr_id() == 14 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 13 );

                // must be 46
                REQUIRE( tElement->get_hmr_id() == 13 );

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

                    // must be 38
                    REQUIRE( tElement->get_hmr_id() == 38 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 1
                    REQUIRE( tElement->get_owner() == 1);

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 13 );

                    // must be 7
                    REQUIRE( tElement->get_hmr_id() == 37 );

                    // must be false
                    REQUIRE( ! tElement->is_active( tActivePattern ) );

                    // must be true
                    REQUIRE( tElement->is_refined( tActivePattern ) );

                    // must be UINT_MAX
                    // (padding elements do not belong to anybody )
                    REQUIRE( tElement->get_owner() == moris::hmr::gNoProcOwner );

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 8 );

                    // must be 32
                    REQUIRE( tElement->get_hmr_id() == 32 );

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
                    moris::Matrix< moris::DDLUMat > tAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 1 );

                    REQUIRE( tAuraElements.length() == 4 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 26 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 27 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 32 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 33 );

                    // ------------------------------------------------------
                    // test aura from proc 0 to proc 1
                    // ------------------------------------------------------

                    // test inverse aura
                    moris::Matrix< moris::DDLUMat > tInverseAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 1 );

                    // must be 4
                    REQUIRE( tInverseAuraElements.length() == 4 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 38 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 39 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 44 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 45 );
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

                    REQUIRE( tElement->get_hmr_id() == 38 );

                    REQUIRE( tElement->is_active( tActivePattern ) );

                    REQUIRE( !tElement->is_refined( tActivePattern ) );

                    REQUIRE( tElement->get_owner() == 2 );

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 8 );

                    // must be 32
                    REQUIRE( tElement->get_hmr_id() == 32 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 1
                    REQUIRE( tElement->get_owner() ==  1);

                    // ------------------------------------------------------

                    // test element
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 32 );

                    REQUIRE( tElement->get_hmr_id() == 56 );

                    // must be true
                    REQUIRE( tElement->is_active( tActivePattern ) );

                    // must be false
                    REQUIRE( ! tElement->is_refined( tActivePattern ) );

                    // must be 3
                    REQUIRE( tElement->get_owner() == 3 );

                    // ------------------------------------------------------
                    // test aura from proc 2 to proc 1
                    // ------------------------------------------------------

                    moris::Matrix< moris::DDLUMat > tAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 1 );

                    // must be 4
                    REQUIRE( tAuraElements.length() == 4);

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 26 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 27 );

                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 32 );
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 33 );

                    // ------------------------------------------------------
                    // test aura from proc 1 to proc 2
                    // ------------------------------------------------------

                    moris::Matrix< moris::DDLUMat > tInverseAuraElements
                    = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 1 );

                    // must be
                    REQUIRE( tInverseAuraElements.length() == 4 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                    REQUIRE( tElement->get_hmr_id() == 38 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                    REQUIRE( tElement->get_hmr_id() == 39 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                    REQUIRE( tElement->get_hmr_id() == 44 );

                    // must be
                    tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                    REQUIRE( tElement->get_hmr_id() == 45 );
                }
            }

            tParameters->set_severity_level( 2 );

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
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 2;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {2}, {8} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 0 );
            tParameters->set_staircase_buffer( 0 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            if ( moris::par_size() == 1 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // since no element has been refined yet, I know that element
                // element 0 is the element with ID 14
                tElement = tBackgroundMesh->get_element( 0 );

                // refine it
                tBackgroundMesh->refine_element( tElement, false );

                // since the list of active elements has not been updated yet,
                // In know that 3 is the element with ID 21
                tElement = tBackgroundMesh->get_element( 3 );

                // refine it
                tBackgroundMesh->refine_element( tElement, false );

                // update element table
                tBackgroundMesh->collect_active_elements();

                //                tBackgroundMesh->save_to_vtk( "BackgorundMesh_refine.vtk");

                // list of active elements on proc
                moris::Matrix< moris::DDLUMat > tElementIDs;

                tBackgroundMesh->get_active_elements_on_proc( tElementIDs );

                REQUIRE( tElementIDs.length()  == 22 );

                REQUIRE( tElementIDs(  0 ) == 124 );
                REQUIRE( tElementIDs(  1 ) == 125 );
                REQUIRE( tElementIDs(  2 ) == 136 );
                REQUIRE( tElementIDs(  3 ) == 137 );
                REQUIRE( tElementIDs(  4 ) ==  15 );
                REQUIRE( tElementIDs(  5 ) ==  20 );
                REQUIRE( tElementIDs(  6 ) == 150 );
                REQUIRE( tElementIDs(  7 ) == 151 );
                REQUIRE( tElementIDs(  8 ) == 162 );
                REQUIRE( tElementIDs(  9 ) == 163 );
                REQUIRE( tElementIDs( 10 ) ==  26 );
                REQUIRE( tElementIDs( 11 ) ==  27 );
                REQUIRE( tElementIDs( 12 ) ==  32 );
                REQUIRE( tElementIDs( 13 ) ==  33 );
                REQUIRE( tElementIDs( 14 ) ==  38 );
                REQUIRE( tElementIDs( 15 ) ==  39 );
                REQUIRE( tElementIDs( 16 ) ==  44 );
                REQUIRE( tElementIDs( 17 ) ==  45 );
                REQUIRE( tElementIDs( 18 ) ==  50 );
                REQUIRE( tElementIDs( 19 ) ==  51 );
                REQUIRE( tElementIDs( 20 ) ==  56 );
                REQUIRE( tElementIDs( 21 ) ==  57 );
            }
            else if ( moris::par_size() == 2 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base * tElement = nullptr;

                // update element table
                tBackgroundMesh->collect_active_elements();

                if( moris::par_rank() == 0 )
                {
                    // pick element 32
                    tElement = tBackgroundMesh->get_element( 6 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 32 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                else
                {
                    // pick element 38
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 38 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement( 0 );

                if( moris::par_rank() == 0 )
                {
                    // pick element 196
                    tElement =  tBackgroundMesh->get_element( 6 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 196 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement( 0 );

                // list of active elements on proc
                moris::Matrix< moris::DDLUMat > tElementIDs;

                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 21 );
                //std::cout<<moris::par_rank()<<": "<<tElementIDs.length()<<std::endl;

                if ( moris::par_rank() == 1 )
                {
                    REQUIRE( tElementIDs(  0 ) ==  26 );
                    REQUIRE( tElementIDs(  1 ) ==  27 );
                    REQUIRE( tElementIDs(  2 ) == 848 );
                    REQUIRE( tElementIDs(  3 ) == 849 );
                    REQUIRE( tElementIDs(  4 ) == 872 );
                    REQUIRE( tElementIDs(  5 ) == 873 );
                    REQUIRE( tElementIDs(  6 ) == 197 );
                    REQUIRE( tElementIDs(  7 ) == 208 );
                    REQUIRE( tElementIDs(  8 ) == 209 );
                    REQUIRE( tElementIDs(  9 ) ==  33 );
                    REQUIRE( tElementIDs( 10 ) == 220 );
                    REQUIRE( tElementIDs( 11 ) == 221 );
                    REQUIRE( tElementIDs( 12 ) == 232 );
                    REQUIRE( tElementIDs( 13 ) == 233 );
                    REQUIRE( tElementIDs( 14 ) ==  39 );
                    REQUIRE( tElementIDs( 15 ) ==  44 );
                    REQUIRE( tElementIDs( 16 ) ==  45 );
                    REQUIRE( tElementIDs( 17 ) ==  50 );
                    REQUIRE( tElementIDs( 18 ) ==  51 );
                    REQUIRE( tElementIDs( 19 ) ==  56 );
                    REQUIRE( tElementIDs( 20 ) ==  57 );
                }

            }
            else if ( moris::par_size() == 4 )
            {
                // create pointer to an element
                moris::hmr::Background_Element_Base * tElement = nullptr;

                // update element table
                tBackgroundMesh->collect_active_elements();

                if( moris::par_rank() == 1 )
                {
                    // pick element 32
                    tElement = tBackgroundMesh->get_element( 2 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 32 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }
                else if ( moris::par_rank()==2 )
                {
                    // pick element 38
                    tElement = tBackgroundMesh->get_element( 0 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 38 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // call refinement procedure
                tBackgroundMesh->perform_refinement( 0 );

                if( moris::par_rank() == 1 )
                {
                    // pick element 209
                    tElement =  tBackgroundMesh->get_element( 5 );

                    // make sure we picked the correct element
                    REQUIRE( tElement->get_hmr_id() == 209 );

                    // flag element for refinement
                    tElement->put_on_refinement_queue();
                }

                // refine all elements
                tBackgroundMesh->perform_refinement( 0 );

                // list of active elements on proc
                moris::Matrix< moris::DDLUMat > tElementIDs;

                if ( moris::par_rank() == 0 )
                {
                    // get list of active elements including aura
                    tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                    // in this case, tElementIDs are different for each proc
                    REQUIRE( tElementIDs.length() == 14 );
                }
                else if ( moris::par_rank() == 1 or moris::par_rank() == 2 )
                {
                    // get list of active elements including aura
                    tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                    // in this case, tElementIDs count is the same for procs 1 and 2
                    REQUIRE( tElementIDs.length() == 21 );
                    if (moris::par_rank()==2)
                    {
                        for (uint i = 0; i < 21; ++i)
                        {
                            REQUIRE( tElementIDs(  0 ) ==  26 );
                            REQUIRE( tElementIDs(  1 ) ==  27 );
                            REQUIRE( tElementIDs(  2 ) == 196 );
                            REQUIRE( tElementIDs(  3 ) == 197 );
                            REQUIRE( tElementIDs(  4 ) == 208 );
                            REQUIRE( tElementIDs(  5 ) == 898 );
                            REQUIRE( tElementIDs(  6 ) == 899 );
                            REQUIRE( tElementIDs(  7 ) == 922 );
                            REQUIRE( tElementIDs(  8 ) == 923 );
                            REQUIRE( tElementIDs(  9 ) ==  33 );
                            REQUIRE( tElementIDs( 10 ) == 220 );
                            REQUIRE( tElementIDs( 11 ) == 221 );
                            REQUIRE( tElementIDs( 12 ) == 232 );
                            REQUIRE( tElementIDs( 13 ) == 233 );
                            REQUIRE( tElementIDs( 14 ) ==  39 );
                            REQUIRE( tElementIDs( 15 ) ==  44 );
                            REQUIRE( tElementIDs( 16 ) ==  45 );
                            REQUIRE( tElementIDs( 17 ) ==  50 );
                            REQUIRE( tElementIDs( 18 ) ==  51 );
                            REQUIRE( tElementIDs( 19 ) ==  56 );
                            REQUIRE( tElementIDs( 20 ) ==  57 );
                        }
                    }
                }
                else if ( moris::par_rank() == 3 )
                {
                    // get list of active elements including aura
                    tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                    // in this case, tElementIDs must be identical on both procs
                    REQUIRE( tElementIDs.length() == 11 );
                }

            }

            // delete background mesh
            delete tBackgroundMesh;

            // delete settings object
            delete tParameters;
        }
    }
}

TEST_CASE("HMR_Background_Mesh_Proc_Nonmultiple_Background", "[moris],[mesh],[hmr],[Background_Nonmultiple],[Background_Mesh]")

{
    if( moris::par_size() == 1 ||  moris::par_size() == 2  || moris::par_size() == 4 )
    {

        // create settings object
        auto tParameters = new moris::hmr::Parameters;

        // Processor Decomposition Method
        uint tDecompMethod = 1;
        tParameters->set_processor_decomp_method(tDecompMethod);

        // set number of elements
        moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {5} };
        tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

        // do not print debug information during test
        tParameters->set_severity_level( 0 );

        // set buffer size to zero
        tParameters->set_refinement_buffer( 1 );
        tParameters->set_staircase_buffer( 1 );

        tParameters->set_lagrange_orders  ( { {2} });
        tParameters->set_lagrange_patterns({ {0} });

        // deactivate truncation
        tParameters->set_bspline_truncation( false );

        // create factory
        moris::hmr::Factory tFactory( tParameters );

        // create background mesh object
        moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

        // update element table
        //tBackgroundMesh->collect_active_elements();

        uint tActivePattern = tBackgroundMesh->get_activation_pattern();

        //            tBackgroundMesh->save_to_vtk( "BackgorundMesh_2Proc.vtk");

        if ( moris::par_size() == 1 )
        {
            // element pointer for tests
            moris::hmr::Background_Element_Base* tElement = nullptr;

            // ------------------------------------------------------
            // simple element test
            // ------------------------------------------------------

            //  tBackgroundMesh->print_level_zero();
            tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 18 );
            // must be 18
            REQUIRE( tElement->get_hmr_id() == 18 );
            // must be true
            REQUIRE( tElement->is_active( tActivePattern ) );
            // must be false
            REQUIRE( ! tElement->is_refined( tActivePattern ) );

            //  tBackgroundMesh->print_level_zero();
            tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 53 );
            // must be 53
            REQUIRE( tElement->get_hmr_id() == 53 );
            // must be true
            REQUIRE( tElement->is_active( tActivePattern ) );
            // must be false
            REQUIRE( ! tElement->is_refined( tActivePattern ) );
        }
        else if (moris::par_size() == 2 )
        {
            if (moris::par_rank() == 0)
            {
                // print zero layer for debugging
                // tBackgroundMesh->print_level_zero();

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 15 );

                // must be 19
                REQUIRE( tElement->get_hmr_id() == 19 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                // must be 1
                REQUIRE( tElement->get_owner() == 0);

                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 16 );

                // must be 20
                REQUIRE( tElement->get_hmr_id() == 20 );

                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );

                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );

                // must be 1
                REQUIRE( tElement->get_owner() == 1);

                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 9 );

                // must be 11
                REQUIRE( tElement->get_hmr_id() == 11 );

                // must be false
                REQUIRE( ! tElement->is_active( tActivePattern ) );

                // must be true
                REQUIRE( tElement->is_refined( tActivePattern ) );

                // must be UINT_MAX
                // (padding elements do not belong to anybody )
                REQUIRE( tElement->get_owner() == moris::hmr::gNoProcOwner );

                // ------------------------------------------------------
                // test aura from proc 1 to proc 0
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 5 );

                REQUIRE( tAuraElements.length() == 10 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 20 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 21 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                REQUIRE( tElement->get_hmr_id() == 28 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                REQUIRE( tElement->get_hmr_id() == 29 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 4 ) );
                REQUIRE( tElement->get_hmr_id() == 36 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 5 ) );
                REQUIRE( tElement->get_hmr_id() == 37 );

                // ------------------------------------------------------
                // test aura from proc 0 to proc 1
                // ------------------------------------------------------

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 5 );

                // must be 10
                REQUIRE( tInverseAuraElements.length() == 10 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 18 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 19 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                REQUIRE( tElement->get_hmr_id() == 26 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                REQUIRE( tElement->get_hmr_id() == 27 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 4 ) );
                REQUIRE( tElement->get_hmr_id() == 34 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 5 ) );
                REQUIRE( tElement->get_hmr_id() == 35 );
            }
        }
        else if (moris::par_size() == 4)
        {
            if (moris::par_rank() == 0)
            {
                // print zero layer for debugging
                // tBackgroundMesh->print_level_zero();

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 27 );
                // must be 35
                REQUIRE( tElement->get_hmr_id() == 35 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 1
                REQUIRE( tElement->get_owner() == 0);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 28 );
                // must be 36
                REQUIRE( tElement->get_hmr_id() == 36 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 2);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 33 );
                // must be 43
                REQUIRE( tElement->get_hmr_id() == 43 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 1);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 34 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 44 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 3);

                // ------------------------------------------------------
                // test aura from proc 0 to proc 1
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 7 );

                REQUIRE( tAuraElements.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 42 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 43 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 2 ) );
                REQUIRE( tElement->get_hmr_id() == 50 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 3 ) );
                REQUIRE( tElement->get_hmr_id() == 51 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 7 );

                REQUIRE( tInverseAuraElements.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 26 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 27 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 2 ) );
                REQUIRE( tElement->get_hmr_id() == 34 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 3 ) );
                REQUIRE( tElement->get_hmr_id() == 35 );

                // ------------------------------------------------------
                // test aura from proc 0 to proc 2
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements2
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 5 );

                REQUIRE( tAuraElements2.length() == 6 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements2( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 20 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements2( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 21 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements2
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 5 );

                REQUIRE( tInverseAuraElements2.length() == 6 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements2( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 18 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements2( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 19 );

                // ------------------------------------------------------
                // test aura from proc 0 to proc 3
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements3
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 8 );

                REQUIRE( tAuraElements3.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements3( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 44 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements3( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 45 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements3
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 8 );

                REQUIRE( tInverseAuraElements3.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements3( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 26 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements3( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 27 );

            }
            else if ( moris::par_rank() == 3)
            {
                // print zero layer for debugging
                // tBackgroundMesh->print_level_zero();

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 7 );
                // must be 35
                REQUIRE( tElement->get_hmr_id() == 35 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 1
                REQUIRE( tElement->get_owner() == 0);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 8 );
                // must be 36
                REQUIRE( tElement->get_hmr_id() == 36 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 2);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 13 );
                // must be 43
                REQUIRE( tElement->get_hmr_id() == 43 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 1);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 14 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 44 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 3);

                // ------------------------------------------------------
                // test aura from proc 3 to proc 0
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 0 );

                REQUIRE( tAuraElements.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 26 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 27 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 0 );

                REQUIRE( tInverseAuraElements.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 44 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 45 );

                // ------------------------------------------------------
                // test aura from proc 3 to proc 1
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements1
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 3 );

                REQUIRE( tAuraElements1.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements1( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 42 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements1( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 43 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements1
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 3 );

                REQUIRE( tInverseAuraElements1.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements1( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 44 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements1( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 45 );

                // ------------------------------------------------------
                // test aura from proc 3 to proc 2
                // ------------------------------------------------------

                // get aura elements
                moris::Matrix< moris::DDLUMat > tAuraElements2
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_aura( 1 );

                REQUIRE( tAuraElements2.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements2( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 28 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tAuraElements2( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 29 );

                // test inverse aura
                moris::Matrix< moris::DDLUMat > tInverseAuraElements2
                = tBackgroundMesh->get_subdomain_ids_of_coarsest_inverse_aura( 1 );

                REQUIRE( tInverseAuraElements2.length() == 4 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements2( 0 ) );
                REQUIRE( tElement->get_hmr_id() == 44 );

                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( tInverseAuraElements2( 1 ) );
                REQUIRE( tElement->get_hmr_id() == 45 );
            }
        }
    }

}

TEST_CASE("HMR_Background_Mesh_Proc_Nonmultiple_Background_3D", "[moris],[mesh],[hmr],[Background_Nonmultiple_3D],[Background_Mesh]")
{
    SECTION( "Background mesh initialization test")
        {
        if( moris::par_size() == 12 )
        {

            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 2;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {5}, {8} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_severity_level( 0 );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 1 );
            tParameters->set_staircase_buffer( 1 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            //tBackgroundMesh->collect_active_elements();

            uint tActivePattern = tBackgroundMesh->get_activation_pattern();

            //tBackgroundMesh->save_to_vtk( "BackgorundMesh_3D.vtk");

            if (moris::par_rank() == 0)
            {
                // print zero layer for debugging
                // tBackgroundMesh->print_level_zero();

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 195 );
                // must be 323
                REQUIRE( tElement->get_hmr_id() == 323 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 1
                REQUIRE( tElement->get_owner() == 0);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 196 );
                // must be 36
                REQUIRE( tElement->get_hmr_id() == 324 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 6);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 201 );
                // must be 43
                REQUIRE( tElement->get_hmr_id() == 331 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 3);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 202 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 332 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 9);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 237 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 395 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 1);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 238 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 396 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 7);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 243 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 403 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 4);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 244 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 404 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 10);

            }
            else if ( moris::par_rank() == 1)
            {
                // print zero layer for debugging
                // tBackgroundMesh->print_level_zero();

                // element pointer for tests
                moris::hmr::Background_Element_Base* tElement = nullptr;

                // ------------------------------------------------------
                // simple element test
                // ------------------------------------------------------

                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 195 );
                // must be 323
                REQUIRE( tElement->get_hmr_id() == 539 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 1
                REQUIRE( tElement->get_owner() == 1);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 196 );
                // must be 36
                REQUIRE( tElement->get_hmr_id() == 540 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 7);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 201 );
                // must be 43
                REQUIRE( tElement->get_hmr_id() == 547 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 4);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 202 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 548 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 10);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 237 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 611 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 2);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 238 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 612 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 8);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 243 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 619 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 5);

                // ------------------------------------------------------
                // test element
                tElement = tBackgroundMesh->get_coarsest_element_by_subdomain_id( 244 );
                // must be 44
                REQUIRE( tElement->get_hmr_id() == 620 );
                // must be true
                REQUIRE( tElement->is_active( tActivePattern ) );
                // must be false
                REQUIRE( ! tElement->is_refined( tActivePattern ) );
                // must be 2
                REQUIRE( tElement->get_owner() == 11);
            }
        }
        }
    SECTION( "Background mesh refinement test")
    {
        if( moris::par_size() == 12 )
        {
            // create settings object
            auto tParameters = new moris::hmr::Parameters;

            // Processor Decomposition Method
            uint tDecompMethod = 2;
            tParameters->set_processor_decomp_method(tDecompMethod);

            // set number of elements
            moris::Matrix< moris::DDLUMat > tNumberOfElements = { {4}, {5}, {8} };
            tParameters->set_number_of_elements_per_dimension( tNumberOfElements );

            // do not print debug information during test
            tParameters->set_severity_level( 0 );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 1 );
            tParameters->set_staircase_buffer( 1 );

            tParameters->set_lagrange_orders  ( { {2} });
            tParameters->set_lagrange_patterns({ {0} });

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // create factory
            moris::hmr::Factory tFactory( tParameters );

            // create background mesh object
            moris::hmr::Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            // update element table
            //tBackgroundMesh->collect_active_elements();

            // create pointer to an element
            moris::hmr::Background_Element_Base * tElement = nullptr;

            // update element table
            tBackgroundMesh->collect_active_elements();

            if( moris::par_rank() == 0 )
            {
                // pick element 17
                tElement = tBackgroundMesh->get_element( 17 );

                // make sure we picked the correct element
                REQUIRE( tElement->get_hmr_id() == 323 );

                // flag element for refinement
                tElement->put_on_refinement_queue();
            }

            // refine all elements
            tBackgroundMesh->perform_refinement( 0 );

            if( moris::par_rank() == 0 )
            {
                // pick refined element 24
                tElement = tBackgroundMesh->get_element( 24 );

                // make sure we picked the correct element
                REQUIRE( tElement->get_hmr_id() == 3607 );

                // flag element for refinement
                tElement->put_on_refinement_queue();
            }

            // refine all elements
            tBackgroundMesh->perform_refinement( 0 );

            //tBackgroundMesh->save_to_vtk( "BackgorundMesh_3D_refine.vtk");

            // list of active elements on proc
            moris::Matrix< moris::DDLUMat > tElementIDs;

            if (moris::par_rank() == 0 or moris::par_rank() == 6)
            {
                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 163 );
            }
            else if (moris::par_rank() == 3 or moris::par_rank() == 9)
            {
                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 143 );
            }
            else if (moris::par_rank() == 1 or moris::par_rank() == 7)
            {
                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 203 );
            }
            else if (moris::par_rank() == 4 or moris::par_rank() == 10)
            {
                // get list of active elements including aura
                tBackgroundMesh->get_active_elements_on_proc_including_aura( tElementIDs );

                // in this case, tElementIDs must be identical on both procs
                REQUIRE( tElementIDs.length() == 175 );
            }

        }
    }
}

