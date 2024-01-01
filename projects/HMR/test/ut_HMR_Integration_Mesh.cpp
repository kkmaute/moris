/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_Integration_Mesh.cpp
 *
 */

#include "catch.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_HMR_Field.hpp"
#include "fn_norm.hpp"

namespace moris::hmr
{
    static real LevelSetFunction1( const Matrix< DDRMat >& aPoint )
    {
        return norm( aPoint ) - 0.5;
    }

    static real LevelSetFunction2( const Matrix< DDRMat >& aPoint )
    {
        return norm( aPoint ) - 4.1;
    }

    TEST_CASE( "HMR Integration Mesh", "[hmr],[IG_Mesh]" )
    {
        //------------------------------------------------------------------------------

        uint tLagrangeMeshIndex = 0;
        uint tBSplineMeshIndex = 0;

        Parameters tParameters;

        tParameters.set_number_of_elements_per_dimension( { { 10 }, { 4 }, { 4 } } );
        tParameters.set_domain_dimensions( { { 10 }, { 4 }, { 4 } } );
        tParameters.set_domain_offset( { { -10.0 }, { -2.0 }, { -2.0 } } );
        tParameters.set_bspline_truncation( true );
        tParameters.set_side_sets( { { 1 }, { 6 }, { 3 }, { 4 }, { 5 }, { 2 } } );

        tParameters.set_output_meshes( { { { 0 } } } );

        tParameters.set_lagrange_orders( { { 1 } } );
        tParameters.set_lagrange_patterns( { { 0 } } );

        tParameters.set_bspline_orders( { { 1 } } );
        tParameters.set_bspline_patterns( { { 0 } } );

        tParameters.set_union_pattern( 2 );
        tParameters.set_working_pattern( 3 );

        tParameters.set_refinement_buffer( 1 );
        tParameters.set_staircase_buffer( 1 );

        Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
        tLagrangeToBSplineMesh( 0 ) = { { 0 } };

        tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

        HMR tHMR( tParameters );

        std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

        // create field
        std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tBSplineMeshIndex );

        // Refine 
        for ( uint iRefinement = 0; iRefinement < 3; iRefinement++ )
        {
            tField->evaluate_scalar_function( LevelSetFunction1 );
            tHMR.flag_surface_elements_on_working_pattern( tField );
            tHMR.perform_refinement_based_on_working_pattern( 0 );
        }

        tHMR.finalize();

        // create pointer to output mesh
        Interpolation_Mesh_HMR* tOutputInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

        // create pointer to output mesh
        Integration_Mesh_HMR* tOutputIntegMesh = tHMR.create_integration_mesh( 1, 0, tOutputInterpMesh );

        Vector< std::string > tBlockNames = tOutputIntegMesh->get_block_set_names();

        Vector< mtk::Cluster const* > tCellClustersInBlock = tOutputIntegMesh->get_cell_clusters_in_set( 0 );

        CHECK( tBlockNames.size() == 1 );
        CHECK( tBlockNames( 0 ) == "HMR_dummy" );
        CHECK( tCellClustersInBlock.size() == tOutputIntegMesh->get_num_elems() );
        CHECK( tOutputInterpMesh->get_num_elems() == tOutputIntegMesh->get_num_elems() );

        uint tSideNames = tOutputIntegMesh->get_num_side_sets();

        CHECK( tSideNames == 6 );

        delete tOutputInterpMesh;
        delete tOutputIntegMesh;
    }

    TEST_CASE( "HMR_Basis_Support", "[hmr][HMR_Basis_Support]" )
    {
        //------------------------------------------------------------------------------

        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );
            tParameters.set_bspline_truncation( true );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 1 } } );

            tParameters.set_refinement_buffer( 1 );
            tParameters.set_staircase_buffer( 1 );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            for ( uint tLevel = 0; tLevel < 1; ++tLevel )
            {
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();

            //        tDatabase->get_background_mesh()->save_to_vtk("Basis_support.vtk");

            tHMR.finalize();

            auto tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            Matrix< IndexMat > tElementIdnInBasisSupport_1;
            Matrix< IndexMat > tElementIdnInBasisSupport_2;
            Matrix< IndexMat > tElementIdnInBasisSupport_3;

            tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 1, tElementIdnInBasisSupport_1 );
            tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 4, tElementIdnInBasisSupport_2 );
            tMesh->get_lagrange_mesh()->get_my_elements_in_basis_support( 0, 8, tElementIdnInBasisSupport_3 );

            CHECK( tElementIdnInBasisSupport_1( 0 ) == 0 );
            CHECK( tElementIdnInBasisSupport_1( 1 ) == 1 );
            CHECK( tElementIdnInBasisSupport_1( 2 ) == 2 );
            CHECK( tElementIdnInBasisSupport_1( 3 ) == 3 );
            CHECK( tElementIdnInBasisSupport_1( 4 ) == 4 );

            CHECK( tElementIdnInBasisSupport_2( 0 ) == 0 );
            CHECK( tElementIdnInBasisSupport_2( 1 ) == 1 );
            CHECK( tElementIdnInBasisSupport_2( 2 ) == 2 );
            CHECK( tElementIdnInBasisSupport_2( 3 ) == 3 );
            CHECK( tElementIdnInBasisSupport_2( 4 ) == 4 );
            CHECK( tElementIdnInBasisSupport_2( 5 ) == 5 );
            CHECK( tElementIdnInBasisSupport_2( 6 ) == 6 );

            CHECK( tElementIdnInBasisSupport_3( 0 ) == 6 );
        }
    }

    TEST_CASE( "HMR Integration Mesh bounding box", "[hmr],[IG_Mesh_bounding_box]" )
    {

        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;
            uint tBSplineMeshIndex = 0;

            Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 10 }, { 4 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 10 }, { 4 }, { 4 } } );
            tParameters.set_domain_offset( { { -5.0 }, { -2.0 }, { -2.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 1 }, { 6 }, { 3 }, { 4 }, { 5 }, { 2 } } );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_refinement_buffer( 1 );
            tParameters.set_staircase_buffer( 1 );

            Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            HMR tHMR( tParameters );

            std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tBSplineMeshIndex );

            // Refine
            for ( uint iRefinement = 0; iRefinement < 2; iRefinement++ )
            {
                tField->evaluate_scalar_function( LevelSetFunction1 );
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );
            }

            tHMR.finalize();

            // create pointer to output mesh
            Interpolation_Mesh_HMR* tOutputInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            Matrix< IndexMat > tNodeIndices;

            tOutputInterpMesh->get_nodes_indices_in_bounding_box( { { 0.1 }, { 1.1 }, { 1.1 } },
                                                                  { { 0.9 }, { 1 }, { 1 } }, tNodeIndices );

            Matrix< IndexMat > tReferenceIndices = { { 241 }, { 301 }, { 430 }, { 431 }, { 639 }, { 638 }, { 747 },
                    { 748 }, { 302 }, { 432 }, { 640 }, { 749 }, { 433 }, { 400 }, { 750 }, { 751 }, { 434 }, { 752 },
                    { 610 }, { 645 }, { 753 }, { 754 }, { 646 }, { 755 }, { 756 }, { 738 }, { 757 }, { 309 }, { 439 },
                    { 649 }, { 758 }, { 310 }, { 440 }, { 650 }, { 759 }, { 441 }, { 760 }, { 442 }, { 761 }, { 653 },
                    { 762 }, { 654 }, { 763 }, { 764 }, { 765 }, { 443 }, { 402 }, { 766 }, { 741 }, { 454 }, { 455 },
                    { 767 }, { 768 }, { 456 }, { 769 }, { 457 }, { 770 }, { 771 }, { 458 }, { 772 }, { 773 }, { 774 },
                    { 775 }, { 776 }, { 777 }, { 616 }, { 657 }, { 778 }, { 744 } };

            REQUIRE( tNodeIndices.numel() >= tReferenceIndices.numel() );
            for ( uint iNode = 0; iNode < tReferenceIndices.numel(); iNode++ )
            {
                CHECK( tNodeIndices( iNode ) == tReferenceIndices( iNode ) );
            }

            delete tOutputInterpMesh;
        }
    }

    TEST_CASE( "HMR delete mesh", "[hmr],[IG_Mesh_delete_mesh]" )
    {
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;
            uint tBSplineMeshIndex = 0;

            Parameters tParameters;

            tParameters.set_number_of_elements_per_dimension( { { 10 }, { 4 }, { 4 } } );
            tParameters.set_domain_dimensions( { { 10 }, { 4 }, { 4 } } );
            tParameters.set_domain_offset( { { -2.0 }, { -2.0 }, { -2.0 } } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_side_sets( { { 1 }, { 6 }, { 3 }, { 4 }, { 5 }, { 2 } } );

            tParameters.set_output_meshes( { { { 0 } } } );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 1 } } );
            tParameters.set_bspline_patterns( { { 0 } } );

            tParameters.set_refinement_buffer( 1 );
            tParameters.set_staircase_buffer( 1 );

            Vector< Matrix< DDSMat >> tLagrangeToBSplineMesh( 1 );
            tLagrangeToBSplineMesh( 0 ) = { { 0 } };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            HMR tHMR( tParameters );

            std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< Field > tField = tMesh->create_field( "Circle", tBSplineMeshIndex );

            tField->evaluate_scalar_function( LevelSetFunction1 );

            // refine
            for ( uint iRefinement = 0; iRefinement < 2; iRefinement++ )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField->evaluate_scalar_function( LevelSetFunction1 );
            }

            // finalize
            tHMR.finalize();

            // tHMR.save_to_exodus( 0, "delete_test_1.exo" );

            // create pointer to output mesh
            Interpolation_Mesh_HMR* tOutputInterpMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            delete tOutputInterpMesh;

            // reset
            tHMR.reset_HMR();

            std::shared_ptr< Mesh > tMesh1 = tHMR.create_mesh( tLagrangeMeshIndex );

            // create field
            std::shared_ptr< Field > tField1 = tMesh1->create_field( "Circle", tBSplineMeshIndex );

            tField1->evaluate_scalar_function( LevelSetFunction2 );

            // refine
            for ( uint iRefinement = 0; iRefinement < 2; iRefinement++ )
            {
                tHMR.flag_surface_elements_on_working_pattern( tField1 );
                tHMR.perform_refinement_based_on_working_pattern( 0 );

                tField1->evaluate_scalar_function( LevelSetFunction2 );
            }

            // finalize
            tHMR.finalize();

            Interpolation_Mesh_HMR* tOutputInterpMesh1 = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // tHMR.save_to_exodus( 0, "delete_test_2.exo" );

            delete tOutputInterpMesh1;
        }
    }

}    // namespace moris
