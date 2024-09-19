/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_BSpline_Mesh.cpp
 *
 */

#include <catch.hpp>
#include "cl_HMR_Background_Mesh.hpp"      //HMR/src
#include "cl_HMR_BSpline_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Element.hpp"              //HMR/src
#include "cl_HMR_Factory.hpp"              //HMR/src
#include "cl_HMR_Parameters.hpp"           //HMR/src

#include "cl_Communication_Tools.hpp"      //COM/src
#include "moris_typedefs.hpp"                    //COR/src
#include "cl_Matrix.hpp"                   //LINALG/src

// This test creates a simple refinement pattern and makes sure that each B-Spline
// is only generated once.

namespace moris::hmr
{

    TEST_CASE( "HMR_Bspline_Mesh", "[moris],[mesh],[hmr],[BsplineMesh]" )
    {
        //-------------------------------------------------------------------------------

        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {
            //-------------------------------------------------------------------------------

            SECTION( "B-Spline Mesh 2D: test basis uniqueness" )
            {
                // create settings object
                auto tParameters = new Parameters;

                // set number of elements
                Vector< luint > tNumberOfElementsPerDimension;
                if ( par_size() == 1 )
                {
                    tNumberOfElementsPerDimension.resize( 2, 3 );
                }
                else if ( par_size() == 2 )
                {
                    tNumberOfElementsPerDimension.resize( 2, 6 );
                }
                else if ( par_size() == 4 )
                {
                    tNumberOfElementsPerDimension.resize( 2, 10 );
                }

                tParameters->set_number_of_elements_per_dimension( tNumberOfElementsPerDimension );

                // deactivate truncation
                tParameters->set_bspline_truncation( false );

                // create factory
                Factory tFactory( tParameters );

                // max level to refine
                uint tMaxLevel = 3;

                // loop over several orders
                for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
                {
                    // set buffer size to zero
                    tParameters->set_refinement_buffer( tOrder );
                    tParameters->set_staircase_buffer( tOrder );

                    // set mesh order
                    tParameters->set_lagrange_orders( { tOrder } );
                    tParameters->set_lagrange_patterns( { 0 } );

                    tParameters->set_bspline_orders( { tOrder } );
                    tParameters->set_bspline_patterns( { 0 } );

                    // create background mesh object
                    Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                    // refine a few elements in the mesh
                    for ( uint l = 0; l < tMaxLevel; ++l )
                    {
                        auto tNumberOfElements = tBackgroundMesh->get_number_of_active_elements_on_proc();

                        // refine every other element
                        for ( luint k = 0; k < tNumberOfElements; k += 3 )
                        {
                            // get element
                            Background_Element_Base* tElement = tBackgroundMesh->get_element( k );

                            // flag element for refinement
                            tElement->put_on_refinement_queue();
                        }

                        // refine mesh
                        tBackgroundMesh->perform_refinement( 0 );
                    }

                    // create B-Spline mesh
                    BSpline_Mesh_Base* tBSplineMesh = tFactory.create_bspline_mesh(
                            tBackgroundMesh,
                            0,
                            tOrder,
                            MORIS_UINT_MAX );

                    // test basis uniqueness
                    REQUIRE( tBSplineMesh->test_for_double_basis() );

                    // tidy up
                    delete tBSplineMesh;
                    delete tBackgroundMesh;
                }

                delete tParameters;
            }
            //-------------------------------------------------------------------------------

            SECTION( "B-Spline Mesh 3D: test basis uniqueness" )
            {
                // create settings object
                auto tParameters = new Parameters;

                // set number of elements
                Vector< luint > tNumberOfElementsPerDimension;

                if ( par_size() == 1 )
                {
                    tNumberOfElementsPerDimension.resize( 3, 3 );
                }
                else if ( par_size() == 2 )
                {
                    tNumberOfElementsPerDimension.resize( 3, 6 );
                }
                else if ( par_size() == 4 )
                {
                    tNumberOfElementsPerDimension.resize( 3, 10 );
                }

                tParameters->set_number_of_elements_per_dimension( tNumberOfElementsPerDimension );

                // deactivate truncation
                tParameters->set_bspline_truncation( false );

                // create factory
                Factory tFactory( tParameters );

                // max level to refine
                uint tMaxLevel = 3;

                // loop over several orders
                for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
                {

                    // set buffer size to zero
                    tParameters->set_refinement_buffer( tOrder );
                    tParameters->set_staircase_buffer( tOrder );

                    tParameters->set_lagrange_orders( { tOrder } );
                    tParameters->set_lagrange_patterns( { 0 } );

                    tParameters->set_bspline_orders( { tOrder } );
                    tParameters->set_bspline_patterns( { 0 } );

                    // create background mesh object
                    Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

                    // refine a few elements in the mesh
                    for ( uint l = 0; l < tMaxLevel; ++l )
                    {
                        auto tNumberOfElements = tBackgroundMesh->get_number_of_active_elements_on_proc();

                        // refine every other element
                        for ( luint k = 0; k < tNumberOfElements; k += 3 )
                        {
                            // get element
                            Background_Element_Base* tElement = tBackgroundMesh->get_element( k );

                            // flag element for refinement
                            tElement->put_on_refinement_queue();
                        }

                        // refine mesh
                        tBackgroundMesh->perform_refinement( 0 );
                    }

                    // create B-Spline mesh
                    BSpline_Mesh_Base* tBSplineMesh = tFactory.create_bspline_mesh( tBackgroundMesh, 0, tOrder, MORIS_UINT_MAX );

                    // test basis uniqueness
                    REQUIRE( tBSplineMesh->test_for_double_basis() );
                    // std::cout << "check " << par_rank() << " " << tOrder << " " << tBSplineMesh->test_for_double_basis() << std::endl;
                    //  tidy up
                    delete tBSplineMesh;
                    delete tBackgroundMesh;
                }

                delete tParameters;
            }
        }
    }

    TEST_CASE( "HMR_Bspline_Mesh_Pattern", "[moris],[mesh],[hmr],[Bspline_mesh_pattern],[BsplineMesh]" )
    {
        //-------------------------------------------------------------------------------

        if ( par_size() == 1 )
        {
            // create settings object
            auto tParameters = new Parameters;

            // set number of elements
            tParameters->set_number_of_elements_per_dimension( { 4, 4 } );

            // deactivate truncation
            tParameters->set_bspline_truncation( false );

            // set buffer size to zero
            tParameters->set_refinement_buffer( 1 );
            tParameters->set_staircase_buffer( 1 );

            tParameters->set_lagrange_orders( { 2 } );
            tParameters->set_lagrange_patterns( { 0 } );

            tParameters->set_bspline_orders( { 2 } );
            tParameters->set_bspline_patterns( { 0 } );

            // create factory
            Factory tFactory( tParameters );

            // create background mesh object
            Background_Mesh_Base* tBackgroundMesh = tFactory.create_background_mesh();

            //----------------------------------------------------------------------------------------------------------
            // Work on activation pattern 0 mesh
            tBackgroundMesh->set_activation_pattern( 0 );

            // element 0 is the element with ID 18
            tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
            tBackgroundMesh->perform_refinement( 0 );

            tBackgroundMesh->get_element( 0 )->put_on_refinement_queue();
            tBackgroundMesh->perform_refinement( 0 );

            //----------------------------------------------------------------------------------------------------------
            // Work on activation pattern 1 mesh
            tBackgroundMesh->set_activation_pattern( 1 );

            // element 0 is the element with ID 18
            tBackgroundMesh->get_element( 15 )->put_on_refinement_queue();
            tBackgroundMesh->perform_refinement( 1 );

            tBackgroundMesh->get_element( 18 )->put_on_refinement_queue();
            tBackgroundMesh->perform_refinement( 1 );

            // create B-Spline mesh
            //( Parameters, Background mesh, Pattern, Order)
            BSpline_Mesh_Base* tBSplineMesh_1 = tFactory.create_bspline_mesh( tBackgroundMesh, 0, 1, MORIS_UINT_MAX );
            BSpline_Mesh_Base* tBSplineMesh_2 = tFactory.create_bspline_mesh( tBackgroundMesh, 1, 1, MORIS_UINT_MAX );

            tBSplineMesh_1->test_sanity();
            tBSplineMesh_2->test_sanity();

            REQUIRE( tBSplineMesh_1->get_number_of_active_basis_on_proc() == 37 );
            REQUIRE( tBSplineMesh_2->get_number_of_active_basis_on_proc() == 37 );

#ifdef MORIS_HAVE_DEBUG
            // Check some basis coordinates of B-Spline mesh 1
            const real* tXYZ_1 = tBSplineMesh_1->get_active_basis( 3 )->get_xyz();
            REQUIRE( tXYZ_1[ 0 ] == 0.0625 );
            REQUIRE( tXYZ_1[ 1 ] == 0.0625 );
            const real* tXYZ_2 = tBSplineMesh_1->get_active_basis( 9 )->get_xyz();
            REQUIRE( tXYZ_2[ 0 ] == 0.25 );
            REQUIRE( tXYZ_2[ 1 ] == 0.125 );
            const real* tXYZ_3 = tBSplineMesh_1->get_active_basis( 27 )->get_xyz();
            REQUIRE( tXYZ_3[ 0 ] == 0.0 );
            REQUIRE( tXYZ_3[ 1 ] == 0.75 );
            const real* tXYZ_4 = tBSplineMesh_1->get_active_basis( 34 )->get_xyz();
            REQUIRE( tXYZ_4[ 0 ] == 0.5 );
            REQUIRE( tXYZ_4[ 1 ] == 1.0 );

            // Check some basis coordinates of B-Spline mesh 2
            const real* tXYZ_5 = tBSplineMesh_2->get_active_basis( 3 )->get_xyz();
            REQUIRE( tXYZ_5[ 0 ] == 0.75 );
            REQUIRE( tXYZ_5[ 1 ] == 0.0 );
            const real* tXYZ_6 = tBSplineMesh_2->get_active_basis( 9 )->get_xyz();
            REQUIRE( tXYZ_6[ 0 ] == 1.0 );
            REQUIRE( tXYZ_6[ 1 ] == 0.25 );
            const real* tXYZ_7 = tBSplineMesh_2->get_active_basis( 27 )->get_xyz();
            REQUIRE( tXYZ_7[ 0 ] == 1.0 );
            REQUIRE( tXYZ_7[ 1 ] == 0.875 );
            const real* tXYZ_8 = tBSplineMesh_2->get_active_basis( 35 )->get_xyz();
            REQUIRE( tXYZ_8[ 0 ] == 0.9375 );
            REQUIRE( tXYZ_8[ 1 ] == 1.0 );
#endif
            delete tBSplineMesh_1;
            delete tBSplineMesh_2;
            delete tBackgroundMesh;
            delete tParameters;
        }
    }
}    // namespace moris::hmr
