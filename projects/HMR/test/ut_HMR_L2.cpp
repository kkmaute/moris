/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_L2.cpp
 *
 */

#include <catch.hpp>
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Mesh_Integration.hpp"
#include "cl_MTK_Field.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"      // COM/src
#include "moris_typedefs.hpp"                    // COR/src
#include "cl_Matrix.hpp"                   // LINALG/src
#include "fn_r2.hpp"
#include "fn_norm.hpp"
#include "fn_equal_to.hpp"

#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::hmr
{

    static real LevelSetFunction( const Matrix< DDRMat >& aPoint )
    {
        return norm( aPoint ) - 1.2;
    }

    TEST_CASE( "HMR_Comm_Table", "[moris],[mesh],[hmr],[hmr_Comm_Table]" )
    {
        if ( par_size() == 4 )
        {
            uint tDimension = 3;
            uint tOrder = 1;

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            Matrix< DDLUMat > tNumberOfElements;

            // set element size
            if ( par_size() == 1 )
            {
                tNumberOfElements.set_size( tDimension, 1, 2 );
            }
            else if ( par_size() == 2 )
            {
                tNumberOfElements.set_size( tDimension, 1, 6 );
            }
            else if ( par_size() == 4 )
            {
                tNumberOfElements.set_size( tDimension, 1, 6 );
            }

            // set values to parameters
            tParameters.set_number_of_elements_per_dimension( tNumberOfElements );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 }, { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

            tParameters.set_bspline_orders( { { 1 }, { 1 } } );
            tParameters.set_bspline_patterns( { { 0 }, { 1 } } );

            tParameters.set_refinement_buffer( 1 );
            tParameters.set_staircase_buffer( 1 );

            tParameters.set_union_pattern( 2 );

            Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
            tLagrangeToBSplineMesh( 0 ) = { 0 };
            tLagrangeToBSplineMesh( 1 ) = { 1 };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            tParameters.set_number_aura( false );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->get_background_mesh()->set_activation_pattern( 0 );

            // refine the first element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                tDatabase->flag_element( 0 );

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            // update database etc
            tDatabase->perform_refinement( 0, false );

            // manually select output pattern
            tDatabase->get_background_mesh()->set_activation_pattern( 1 );

            // refine the last element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc() - 1 );

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }
            // update database etc
            tDatabase->perform_refinement( 1, false );

            // manually create union
            tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

            // update background mesh
            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tOrder );

            Matrix< IdMat > tCommTable = tMesh->get_communication_table();

            if ( par_rank() == 0 )
            {
                CHECK( equal_to( tCommTable( 0 ), 1 ) );
                CHECK( equal_to( tCommTable( 1 ), 2 ) );
                CHECK( equal_to( tCommTable( 2 ), 3 ) );
            }
            if ( par_rank() == 1 )
            {
                CHECK( equal_to( tCommTable( 0 ), 0 ) );
                CHECK( equal_to( tCommTable( 1 ), 3 ) );
            }
            if ( par_rank() == 2 )
            {
                CHECK( equal_to( tCommTable( 0 ), 0 ) );
                CHECK( equal_to( tCommTable( 1 ), 3 ) );
            }
            if ( par_rank() == 3 )
            {
                CHECK( equal_to( tCommTable( 0 ), 0 ) );
                CHECK( equal_to( tCommTable( 1 ), 1 ) );
                CHECK( equal_to( tCommTable( 2 ), 2 ) );
            }
        }
    }

    TEST_CASE( "HMR_CommTable2", "[moris],[mesh],[hmr],[HMR_CommTable2]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        if ( par_size() == 4 )
        {

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameter_List tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "6, 6, 6" );
            tParameters.set( "domain_dimensions", "6, 6, 6" );
            tParameters.set( "domain_offset", "-3.0, -3.0, -3.0" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "union_pattern", 2 );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "use_multigrid", 0 );

            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // tDatabase->get_background_mesh()->save_to_vtk("Background555.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            // print(tDatabase->get_proc_neighbors(), "ProcNeighbors");
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern", "[moris],[mesh],[hmr],[hmr_L2_pattern]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        if ( par_size() == 1 )
        {

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameter_List tParameters = prm::create_hmr_parameter_list();

            tParameters.set( "number_of_elements_per_dimension", "2, 2" );

            tParameters.set( "truncate_bsplines", 1 );
            tParameters.set( "lagrange_orders", "1, 1" );
            tParameters.set( "lagrange_pattern", "0, 1" );
            tParameters.set( "bspline_orders", "1, 1, 1, 2" );
            tParameters.set( "bspline_pattern", "0, 0, 1, 1" );

            tParameters.set( "union_pattern", 2 );

            tParameters.set( "lagrange_to_bspline", "0, 1; 2, 3" );

            tParameters.set( "use_multigrid", 0 );

            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            // refine the first element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                //            tDatabase->flag_element( 0 );
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            tDatabase->get_background_mesh()->save_to_vtk( "Background444.vtk" );

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element three times
            for ( uint tLevel = 0; tLevel < 2; ++tLevel )
            {
                //            tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc()-1 );
                tDatabase->get_background_mesh()->get_element(
                        tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            // manually create union
            tDatabase->unite_patterns( 0, 1, 2 );

            tDatabase->get_background_mesh()->save_to_vtk( "Background555.vtk" );
            // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            //------------------------------------------------------------------------------
            //  Fields
            //------------------------------------------------------------------------------

            // create pointer to input field
            auto tField = tHMR.create_field( "LevelSet", 0, 0 );

            // evaluate function at nodes
            tField->evaluate_scalar_function( LevelSetFunction );

            uint tOutputOrder = 1;

            // map input to output
            tHMR.map_field_to_output_union( tField, tOutputOrder );

            auto tOutputMesh = tHMR.create_mesh( tOutputOrder, 2 );    // order , pattern

            // calculate exact value
            auto tExact = tOutputMesh->create_field( "Exact", 0 );

            tExact->evaluate_scalar_function( LevelSetFunction );
            //------------------------------------------------------------------------------
            //   Test error
            //------------------------------------------------------------------------------

            // determine coefficient of determination
            real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

            // perform test
            REQUIRE( tR2 > 0.99 );
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern3", "[moris],[mesh],[hmr],[hmr_L2_pattern3]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        if ( par_size() == 1 )
        {

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            Matrix< DDLUMat > tNumberOfElements;

            // set values to parameters
            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 2 }, { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

            tParameters.set_bspline_orders( { { 1 }, { 1 }, { 1 }, { 2 } } );
            tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

            tParameters.set_union_pattern( 2 );

            Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
            tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
            tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            tParameters.set_number_aura( false );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            // refine the first element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                //            tDatabase->flag_element( 0 );
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            tDatabase->get_background_mesh()->save_to_vtk( "Background444.vtk" );

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element three times
            for ( uint tLevel = 0; tLevel < 2; ++tLevel )
            {
                //            tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc()-1 );
                tDatabase->get_background_mesh()->get_element(
                        tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            // manually create union
            tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

            tDatabase->get_background_mesh()->save_to_vtk( "Background555.vtk" );
            // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            //------------------------------------------------------------------------------
            //  Fields
            //------------------------------------------------------------------------------

            // create pointer to input field
            auto tField = tHMR.create_field( "LevelSet", 0, 0 );

            // evaluate function at nodes
            tField->evaluate_scalar_function( LevelSetFunction );

            uint tOutputOrder = 2;

            // map input to output
            tHMR.map_field_to_output_union( tField, tOutputOrder );

            auto tOutputMesh = tHMR.create_mesh( tOutputOrder, tParameters.get_union_pattern() );    // order , pattern

            // calculate exact value
            auto tExact = tOutputMesh->create_field( "Exact", 0 );

            tExact->evaluate_scalar_function( LevelSetFunction );
            //------------------------------------------------------------------------------
            //   Test error
            //------------------------------------------------------------------------------

            // determine coefficient of determination
            real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

            // perform test
            REQUIRE( tR2 > 0.99 );
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern4", "[moris],[mesh],[hmr],[hmr_L2_pattern4]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        if ( par_size() == 1 )
        {

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            Matrix< DDLUMat > tNumberOfElements;

            // set values to parameters
            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 }, { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

            tParameters.set_bspline_orders( { { 1 }, { 1 }, { 1 }, { 2 } } );
            tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

            tParameters.set_union_pattern( 2 );

            tParameters.set_staircase_buffer( 2 );

            Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
            tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
            tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            tParameters.set_number_aura( false );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            // refine the first element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                //            tDatabase->flag_element( 0 );
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            tDatabase->get_background_mesh()->save_to_vtk( "Background444.vtk" );

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element two times
            for ( uint tLevel = 0; tLevel < 2; ++tLevel )
            {
                //            tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc()-1 );
                tDatabase->get_background_mesh()->get_element(
                        tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            // manually create union
            tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

            tDatabase->get_background_mesh()->save_to_vtk( "Background555.vtk" );
            // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            //------------------------------------------------------------------------------
            //  Fields
            //------------------------------------------------------------------------------

            // create pointer to input field
            auto tField = tHMR.create_field( "LevelSet", 0, 0 );

            // evaluate function at nodes
            tField->evaluate_scalar_function( LevelSetFunction );

            uint tOutputOrder = 2;

            // map input to output
            tHMR.map_field_to_output_union( tField, tOutputOrder );

            auto tOutputMesh = tHMR.create_mesh( tOutputOrder, tParameters.get_union_pattern() );    // order , pattern

            // calculate exact value
            auto tExact = tOutputMesh->create_field( "Exact", 0 );

            tExact->evaluate_scalar_function( LevelSetFunction );
            //-------------------------------------------------------------------------------------

            //         auto tOutputMesh = tHMR.create_mesh( 1 );        // mesh index
            //
            //         // calculate exact value
            //         auto tExact = tOutputMesh->create_field( "Exact", 1 );
            //
            //         tExact->evaluate_scalar_function( LevelSetFunction );
            //------------------------------------------------------------------------------
            //   Test error
            //------------------------------------------------------------------------------

            // determine coefficient of determination
            real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

            // perform test
            REQUIRE( tR2 > 0.99 );
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern2", "[moris],[mesh],[hmr],[hmr_L2_pattern2]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        if ( par_size() == 1 )
        {
            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            Matrix< DDLUMat > tNumberOfElements;

            // set values to parameters
            tParameters.set_number_of_elements_per_dimension( { { 2 }, { 2 } } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 2 }, { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

            tParameters.set_bspline_orders( { { 1 }, { 1 }, { 1 }, { 2 } } );
            tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

            tParameters.set_union_pattern( 2 );

            tParameters.set_staircase_buffer( 2 );

            Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
            tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
            tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

            tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

            tParameters.set_number_aura( false );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            // std::shared_ptr< Database >
            auto tDatabase = tHMR.get_database();

            // manually select output pattern
            tDatabase->set_activation_pattern( 0 );

            // refine the first element three times
            for ( uint tLevel = 0; tLevel < 3; ++tLevel )
            {
                //            tDatabase->flag_element( 0 );
                tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 0 );
            }

            tDatabase->get_background_mesh()->save_to_vtk( "Background444.vtk" );

            // manually select output pattern
            tDatabase->set_activation_pattern( 1 );

            // refine the last element two times
            for ( uint tLevel = 0; tLevel < 2; ++tLevel )
            {
                //            tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc()-1 );
                tDatabase->get_background_mesh()->get_element(
                        tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                // manually refine, do not reset pattern
                tDatabase->get_background_mesh()->perform_refinement( 1 );
            }

            // manually create union
            tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

            tDatabase->get_background_mesh()->save_to_vtk( "Background555.vtk" );
            // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

            tDatabase->update_bspline_meshes();
            tDatabase->update_lagrange_meshes();
            // calculate T-Matrices etc
            tDatabase->finalize();

            //------------------------------------------------------------------------------
            //  Fields
            //------------------------------------------------------------------------------

            // create pointer to input field
            auto tField = tHMR.create_field( "LevelSet", 0, 0 );

            // evaluate function at nodes
            tField->evaluate_scalar_function( LevelSetFunction );

            uint tOutputOrder = 1;

            // map input to output
            tHMR.map_field_to_output_union( tField, tOutputOrder );

            auto tOutputMesh = tHMR.create_mesh( tOutputOrder, tParameters.get_union_pattern() );    // order , pattern

            // calculate exact value
            auto tExact = tOutputMesh->create_field( "Exact", 0 );

            tExact->evaluate_scalar_function( LevelSetFunction );
            //------------------------------------------------------------------------------
            //   Test error
            //------------------------------------------------------------------------------

            // determine coefficient of determination
            real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

            // perform test
            REQUIRE( tR2 > 0.99 );
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern5", "[moris],[mesh],[hmr],[hmr_L2_pattern5]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {
            for ( uint tDimension = 2; tDimension <= 3; ++tDimension )
            {
                for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
                {
                    //------------------------------------------------------------------------------
                    //  HMR Parameters setup
                    //------------------------------------------------------------------------------

                    // The parameter object controls the behavior of HMR.
                    Parameters tParameters;

                    // set values to parameters
                    Matrix< DDLUMat > tNumberOfElements;

                    // set element size
                    if ( par_size() == 1 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 2 );
                    }
                    else if ( par_size() == 2 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }
                    else if ( par_size() == 4 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }

                    tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

                    // B-Spline truncation is turned on by default.
                    // It is recommended to leave this setting as is.
                    tParameters.set_bspline_truncation( true );

                    tParameters.set_lagrange_orders( { { tOrder }, { tOrder } } );
                    tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

                    tParameters.set_bspline_orders( { { tOrder }, { 2 }, { tOrder }, { 2 } } );
                    tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

                    tParameters.set_union_pattern( 2 );

                    tParameters.set_staircase_buffer( tOrder );

                    Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
                    tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
                    tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

                    tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                    tParameters.set_number_aura( false );

                    //------------------------------------------------------------------------------
                    //  HMR Initialization
                    //------------------------------------------------------------------------------

                    // create the HMR object by passing the settings to the constructor
                    HMR tHMR( tParameters );

                    // std::shared_ptr< Database >
                    auto tDatabase = tHMR.get_database();

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 0 );

                    // refine the first element three times
                    for ( uint tLevel = 0; tLevel < 3; ++tLevel )
                    {
                        //            tDatabase->flag_element( 0 );
                        tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 0 );
                    }

                    //           tDatabase->get_background_mesh()->save_to_vtk("Background444.vtk");

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 1 );

                    // refine the last element two times
                    for ( uint tLevel = 0; tLevel < 2; ++tLevel )
                    {
                        tDatabase->get_background_mesh()->get_element(
                                tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 1 );
                    }

                    // manually create union
                    // tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

                    // tDatabase->get_background_mesh()->save_to_vtk("Background555.vtk");
                    // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

                    tDatabase->update_bspline_meshes();
                    tDatabase->update_lagrange_meshes();
                    // calculate T-Matrices etc
                    tDatabase->finalize();

                    //------------------------------------------------------------------------------
                    //  Fields
                    //------------------------------------------------------------------------------

                    // create pointer to input field
                    auto tField = tHMR.create_field( "LevelSet", 0, 0 );

                    // evaluate function at nodes
                    tField->evaluate_scalar_function( LevelSetFunction );

                    uint tOutputMeshIndex = 1;
                    uint tBsplineMeshIndex = 0;
                    // map input to output
                    tHMR.map_field_to_output( tField, tOutputMeshIndex, tBsplineMeshIndex );

                    // tHMR.save_to_exodus( tOutputMeshIndex, "LevelSet2111.exo" );

                    auto tOutputMesh = tHMR.create_mesh( 1 );    // Mesh index

                    // calculate exact value
                    auto tExact = tOutputMesh->create_field( "Exact", 0 );

                    tExact->evaluate_scalar_function( LevelSetFunction );
                    //------------------------------------------------------------------------------
                    //   Test error
                    //------------------------------------------------------------------------------

                    // determine coefficient of determination
                    real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

                    // perform test
                    REQUIRE( tR2 > 0.99 );
                }
            }
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern6", "[moris],[mesh],[hmr],[hmr_L2_pattern6]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {
            for ( uint tDimension = 2; tDimension <= 3; ++tDimension )
            {
                for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
                {
                    //------------------------------------------------------------------------------
                    //  HMR Parameters setup
                    //------------------------------------------------------------------------------

                    // The parameter object controls the behavior of HMR.
                    Parameters tParameters;

                    // set values to parameters
                    Matrix< DDLUMat > tNumberOfElements;

                    // set element size
                    if ( par_size() == 1 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 2 );
                    }
                    else if ( par_size() == 2 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }
                    else if ( par_size() == 4 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }

                    tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

                    // B-Spline truncation is turned on by default.
                    // It is recommended to leave this setting as is.
                    tParameters.set_bspline_truncation( true );

                    tParameters.set_lagrange_orders( { { tOrder }, { tOrder } } );
                    tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

                    tParameters.set_bspline_orders( { { 1 }, { tOrder }, { 1 }, { tOrder } } );
                    tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

                    tParameters.set_union_pattern( 2 );

                    tParameters.set_staircase_buffer( tOrder );

                    Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
                    tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
                    tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

                    tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                    tParameters.set_number_aura( false );

                    //------------------------------------------------------------------------------
                    //  HMR Initialization
                    //------------------------------------------------------------------------------

                    // create the HMR object by passing the settings to the constructor
                    HMR tHMR( tParameters );

                    // std::shared_ptr< Database >
                    auto tDatabase = tHMR.get_database();

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 0 );

                    // refine the first element three times
                    for ( uint tLevel = 0; tLevel < 3; ++tLevel )
                    {
                        //            tDatabase->flag_element( 0 );
                        tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 0 );
                    }

                    tDatabase->get_background_mesh()->save_to_vtk( "Background444.vtk" );

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 1 );

                    // refine the last element two times
                    for ( uint tLevel = 0; tLevel < 2; ++tLevel )
                    {
                        tDatabase->get_background_mesh()->get_element(
                                tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 1 );
                    }

                    // manually create union
                    // tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

                    // tDatabase->get_background_mesh()->save_to_vtk("Background555.vtk");
                    // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

                    tDatabase->update_bspline_meshes();
                    tDatabase->update_lagrange_meshes();
                    // calculate T-Matrices etc
                    tDatabase->finalize();

                    //------------------------------------------------------------------------------
                    //  Fields
                    //------------------------------------------------------------------------------

                    // create pointer to input field
                    auto tField = tHMR.create_field( "LevelSet", 0, 1 );

                    // evaluate function at nodes
                    tField->evaluate_scalar_function( LevelSetFunction );

                    uint tOutputMeshIndex = 1;
                    uint tBsplineMeshIndex = 1;

                    // map input to output
                    tHMR.map_field_to_output( tField, tOutputMeshIndex, tBsplineMeshIndex );

                    // tHMR.save_to_exodus( tOutputMeshIndex, "LevelSet2111.exo" );

                    auto tOutputMesh = tHMR.create_mesh( 1 );    // mesh index

                    // calculate exact value
                    auto tExact = tOutputMesh->create_field( "Exact", 1 );

                    tExact->evaluate_scalar_function( LevelSetFunction );
                    //------------------------------------------------------------------------------
                    //   Test error
                    //------------------------------------------------------------------------------

                    // determine coefficient of determination
                    real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

                    // perform test
                    REQUIRE( tR2 > 0.99 );
                }
            }
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern7", "[moris],[mesh],[hmr],[hmr_L2_pattern7]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {
            for ( uint tDimension = 2; tDimension <= 3; ++tDimension )
            {
                for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
                {
                    //------------------------------------------------------------------------------
                    //  HMR Parameters setup
                    //------------------------------------------------------------------------------

                    // The parameter object controls the behavior of HMR.
                    Parameters tParameters;

                    // set values to parameters
                    Matrix< DDLUMat > tNumberOfElements;

                    // set element size
                    if ( par_size() == 1 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 2 );
                    }
                    else if ( par_size() == 2 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }
                    else if ( par_size() == 4 )
                    {
                        tNumberOfElements.set_size( tDimension, 1, 6 );
                    }

                    tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

                    // B-Spline truncation is turned on by default.
                    // It is recommended to leave this setting as is.
                    tParameters.set_bspline_truncation( true );

                    tParameters.set_lagrange_orders( { { 1 }, { tOrder } } );
                    tParameters.set_lagrange_patterns( { { 0 }, { 1 } } );

                    tParameters.set_bspline_orders( { { 1 }, { 1 }, { 1 }, { tOrder } } );
                    tParameters.set_bspline_patterns( { { 0 }, { 0 }, { 1 }, { 1 } } );

                    tParameters.set_union_pattern( 2 );

                    tParameters.set_staircase_buffer( tOrder );

                    Vector< Vector< uint > > tLagrangeToBSplineMesh( 2 );
                    tLagrangeToBSplineMesh( 0 ) = { 0, 1 };
                    tLagrangeToBSplineMesh( 1 ) = { 2, 3 };

                    tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                    tParameters.set_number_aura( false );

                    //------------------------------------------------------------------------------
                    //  HMR Initialization
                    //------------------------------------------------------------------------------

                    // create the HMR object by passing the settings to the constructor
                    HMR tHMR( tParameters );

                    // std::shared_ptr< Database >
                    auto tDatabase = tHMR.get_database();

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 0 );

                    // refine the first element three times
                    for ( uint tLevel = 0; tLevel < 3; ++tLevel )
                    {
                        //            tDatabase->flag_element( 0 );
                        tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 0 );
                    }

                    // tDatabase->get_background_mesh()->save_to_vtk("BackgroundPresi1.vtk");

                    // manually select output pattern
                    tDatabase->set_activation_pattern( 1 );

                    // refine the last element two times
                    for ( uint tLevel = 0; tLevel < 2; ++tLevel )
                    {
                        //            tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc()-1 );
                        tDatabase->get_background_mesh()->get_element(
                                tDatabase->get_number_of_elements_on_proc() - 1 )->put_on_refinement_queue();

                        // manually refine, do not reset pattern
                        tDatabase->get_background_mesh()->perform_refinement( 1 );
                    }

                    // manually create union
                    // tDatabase->unite_patterns( 0, 1, tParameters.get_union_pattern() );

                    // tDatabase->get_background_mesh()->save_to_vtk("BackgroundPresi2.vtk");
                    // tHMR.mBSplineMeshes( 1 )->save_to_vtk("BSpline.vtk");

                    tDatabase->update_bspline_meshes();
                    tDatabase->update_lagrange_meshes();
                    // calculate T-Matrices etc
                    tDatabase->finalize();

                    //------------------------------------------------------------------------------
                    //  Fields
                    //------------------------------------------------------------------------------

                    // create pointer to input field
                    auto tField = tHMR.create_field( "LevelSet", 0, 1 );

                    // evaluate function at nodes
                    tField->evaluate_scalar_function( LevelSetFunction );

                    uint tOutputMeshIndex = 1;
                    uint tBsplineMeshIndex = 1;

                    // map input to output
                    tHMR.map_field_to_output( tField, tOutputMeshIndex, tBsplineMeshIndex );

                    // tHMR.save_to_exodus( tOutputMeshIndex, "LevelSetPresi.exo" );

                    auto tOutputMesh = tHMR.create_mesh( 1 );    // mesh index

                    // calculate exact value
                    auto tExact = tOutputMesh->create_field( "Exact", 1 );

                    tExact->evaluate_scalar_function( LevelSetFunction );
                    //------------------------------------------------------------------------------
                    //   Test error
                    //------------------------------------------------------------------------------

                    // determine coefficient of determination
                    real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

                    // perform test
                    if ( par_size() == 1 )
                    {
                        REQUIRE( tR2 > 0.95 );
                    }
                    else if ( par_size() == 2 )
                    {
                        REQUIRE( tR2 > 0.96 );
                    }
                    else if ( par_size() == 4 )
                    {
                        REQUIRE( tR2 > 0.94 );
                    }
                }
            }
        }
    }

    TEST_CASE( "HMR_L2_Test_Pattern8", "[moris],[mesh],[hmr],[hmr_L2_pattern8]" )
    {
        //    gLogger.set_severity_level( 0 );
        // can only perform test for 1, 2 or 4 procs
        // do this test for 2 and 3 dimensions
        if ( par_size() == 1 )
        {
            uint tDimension = 2;

            //------------------------------------------------------------------------------
            //  HMR Parameters setup
            //------------------------------------------------------------------------------

            // The parameter object controls the behavior of HMR.
            Parameters tParameters;

            // set values to parameters
            Matrix< DDLUMat > tNumberOfElements;

            // set element size
            tNumberOfElements.set_size( tDimension, 1, 2 );

            tParameters.set_number_of_elements_per_dimension( { tNumberOfElements } );

            // B-Spline truncation is turned on by default.
            // It is recommended to leave this setting as is.
            tParameters.set_bspline_truncation( true );

            tParameters.set_lagrange_orders( { { 1 } } );
            tParameters.set_lagrange_patterns( { { 0 } } );

            tParameters.set_bspline_orders( { { 2 } } );
            tParameters.set_bspline_patterns( { { 1 } } );

            tParameters.set_union_pattern( 2 );

            tParameters.set_staircase_buffer( 2 );

            tParameters.set_initial_refinement( { { 1 } } );
            tParameters.set_initial_refinement_patterns( { { 0 } } );

            tParameters.set_number_aura( false );

            //------------------------------------------------------------------------------
            //  HMR Initialization
            //------------------------------------------------------------------------------

            // create the HMR object by passing the settings to the constructor
            HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            auto tDatabase = tHMR.get_database();

            // calculate T-Matrices etc
            tDatabase->finalize();

            //------------------------------------------------------------------------------
            //  Fields
            //------------------------------------------------------------------------------

            // create pointer to input field
            auto tField = tHMR.create_field( "LevelSet", 0, 0 );

            // evaluate function at nodes
            tField->evaluate_scalar_function( LevelSetFunction );

            Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( 0 );

            // Create integration mesh
            mtk::Integration_Mesh* tIntegrationMesh = tHMR.create_integration_mesh( 0, tInterpolationMesh );

            mtk::Mesh_Pair tMeshPair( tInterpolationMesh, tIntegrationMesh );

            mtk::Field_Discrete tField_proxy( tMeshPair, 0 );

            tField_proxy.unlock_field();
            tField_proxy.set_values( tField->get_node_values() );

            mtk::Mapper tMapper;

            tMapper.perform_mapping( &tField_proxy, mtk::EntityRank::NODE, mtk::EntityRank::BSPLINE );

            tField->get_coefficients() = tField_proxy.get_coefficients();

            tField->evaluate_nodal_values( tField->get_coefficients() );

            tHMR.save_to_exodus( 0, "LevelSetPresi.exo" );
        }
    }
}
