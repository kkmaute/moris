/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_IO.cpp
 *
 */

#include <catch.hpp>
#include "paths.hpp"

#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_BSpline_Mesh_Base.hpp"       //HMR/src
#include "cl_HMR_Factory.hpp"                 //HMR/src
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"    //HMR/src
#include "cl_HMR_Parameters.hpp"            //HMR/src

#include "cl_Communication_Manager.hpp"    //COM/src
#include "cl_Communication_Tools.hpp"      //COM/src

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_r2.hpp"

namespace moris::hmr
{

    static real LevelSetFunction( const Matrix< DDRMat >& aPoint )
    {
        return norm( aPoint ) - 1.2;
    }

    TEST_CASE( "HMR_IO", "[moris],[hmr],[HMR_IO]" )
    {
        if ( par_size() == 1 )
        {
            for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
            {
                // create settings object
                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 } } );

                tParameters.set_domain_dimensions( { { 1 }, { 1 } } );
                tParameters.set_domain_offset( { { -0.5 }, { -0.5 } } );

                // B-Spline truncation is turned on by default.
                // It is recommended to leave this setting as is.
                tParameters.set_bspline_truncation( true );

                tParameters.set_lagrange_orders( { tOrder } );
                tParameters.set_lagrange_patterns( { 2 } );

                tParameters.set_bspline_orders( { tOrder, tOrder } );
                tParameters.set_bspline_patterns( { 0, 1 } );

                tParameters.set_staircase_buffer( 3 );
                tParameters.set_refinement_buffer( 3 );

                tParameters.set_initial_refinement( { { 1 } } );
                tParameters.set_initial_refinement_patterns( { { 0 } } );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 }, { 1 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------
                //  HMR Initialization
                //------------------------------------------------------------------------------

                // create the HMR object by passing the settings to the constructor
                HMR tHMR( tParameters );

                // std::shared_ptr< Database >
                auto tDatabase = tHMR.get_database();

                // manually select output pattern
                tDatabase->set_activation_pattern( 0 );

                tHMR.perform_initial_refinement();

                // manually select output pattern
                tDatabase->set_activation_pattern( 1 );

                // refine the last element three times
                for ( uint tLevel = 0; tLevel < 3; ++tLevel )
                {
                    tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                    // manually refine, do not reset pattern
                    tDatabase->get_background_mesh()->perform_refinement( 1 );
                }

                tDatabase->unite_patterns( 0, 1, 2 );

                tDatabase->update_bspline_meshes();
                tDatabase->update_lagrange_meshes();
                // calculate T-Matrices etc
                tDatabase->finalize();

                tHMR.save_to_hdf5( "Mesh_Data_test.hdf5", 0 );

                HMR tHMR_Input( "Mesh_Data_test.hdf5" );

                // empty container for B-Spline meshes
                Vector< BSpline_Mesh_Base* > tBSplineMeshes;

                // create factory
                Factory tFactory( tHMR_Input.get_database()->get_parameters() );

                // create first order Lagrange mesh
                Lagrange_Mesh_Base* tLagrangeMesh = tFactory.create_lagrange_mesh(
                        tHMR_Input.get_database()->get_background_mesh(), tBSplineMeshes, 2, 1 );

                REQUIRE( tLagrangeMesh->get_number_of_nodes_on_proc() == 116 );

                // Check some basis coordinates of Lagrange mesh 1
                const real* tXYZ_1 = tLagrangeMesh->get_node_by_index( 2 )->get_xyz();
                REQUIRE( tXYZ_1[ 0 ] == -0.46875 );
                REQUIRE( tXYZ_1[ 1 ] == -0.46875 );
                const real* tXYZ_2 = tLagrangeMesh->get_node_by_index( 71 )->get_xyz();
                REQUIRE( tXYZ_2[ 0 ] == 0.0 );
                REQUIRE( tXYZ_2[ 1 ] == 0.0 );
                const real* tXYZ_3 = tLagrangeMesh->get_node_by_index( 109 )->get_xyz();
                REQUIRE( tXYZ_3[ 0 ] == 0.25 );
                REQUIRE( tXYZ_3[ 1 ] == 0.375 );

                delete tLagrangeMesh;
            }
        }
    }

    TEST_CASE( "HMR_Field_IO_EXO", "[moris],[hmr],[HMR_Field_IO_Exo]" )
    {
        if ( par_size() == 1 )
        {
            for ( uint tOrder = 1; tOrder <= 1; tOrder++ )
            {
                std::string tPrefix = get_base_moris_dir();
                std::string tMeshFileName = tPrefix + "projects/HMR/test/data/hmr_Input_Mesh_Data_1.hdf5";

                HMR tHMR( tMeshFileName );

                std::string tFieldFileName = tPrefix + "projects/HMR/test/data/hmr_Field_IO_Test.exo";

                std::shared_ptr< Field > tField = tHMR.load_field_from_file( "Circle", tFieldFileName, 0, 0 );

                //------------------------------------------------------------------------------
                //   Test error
                //------------------------------------------------------------------------------

                auto tOutputMesh = tHMR.create_mesh( 0 );

                // calculate exact value
                auto tExact = tOutputMesh->create_field( "Exact", 0 );

                tExact->evaluate_scalar_function( LevelSetFunction );

                // determine coefficient of determination
                real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

                // perform test
                REQUIRE( tR2 > 0.99 );
            }
        }
    }

    TEST_CASE( "HMR_Field_IO_HDF", "[moris],[hmr],[HMR_Field_IO_HDF]" )
    {
        if ( par_size() == 1 )
        {
            for ( uint tOrder = 1; tOrder <= 3; tOrder++ )
            {
                uint tLagrangeMeshIndex = 0;
                uint tLBSplineMeshIndex = 2;

                // create settings object
                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 4 } } );

                tParameters.set_domain_dimensions( { { 1 }, { 1 } } );
                tParameters.set_domain_offset( { { -0.5 }, { -0.5 } } );

                // B-Spline truncation is turned on by default.
                // It is recommended to leave this setting as is.
                tParameters.set_bspline_truncation( true );

                tParameters.set_lagrange_orders( { tOrder } );
                tParameters.set_lagrange_patterns( { 2 } );

                tParameters.set_bspline_orders( { 1, 1, 1 } );
                tParameters.set_bspline_patterns( { 0, 1, 2 } );

                tParameters.set_output_meshes( { { { 0 } } } );
                tParameters.set_lagrange_input_mesh( { { 0 } } );

                tParameters.set_staircase_buffer( 3 );
                tParameters.set_refinement_buffer( 3 );

                tParameters.set_initial_refinement( { { 1 } } );
                tParameters.set_initial_refinement_patterns( { { 0 } } );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 }, { 1 }, { 2 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------
                //  HMR Initialization
                //------------------------------------------------------------------------------

                // create the HMR object by passing the settings to the constructor
                HMR tHMR( tParameters );

                // std::shared_ptr< Database >
                auto tDatabase = tHMR.get_database();

                // manually select output pattern
                tDatabase->set_activation_pattern( 0 );

                tHMR.perform_initial_refinement();

                // manually select output pattern
                tDatabase->set_activation_pattern( 1 );

                // refine the last element three times
                for ( uint tLevel = 0; tLevel < 3; ++tLevel )
                {
                    tDatabase->get_background_mesh()->get_element( 0 )->put_on_refinement_queue();

                    // manually refine, do not reset pattern
                    tDatabase->get_background_mesh()->perform_refinement( 1 );
                }

                tDatabase->unite_patterns( 0, 1, 2 );

                tDatabase->update_bspline_meshes();
                tDatabase->update_lagrange_meshes();
                // calculate T-Matrices etc
                tDatabase->finalize();

                std::string tPrefix = get_base_moris_dir();

                std::string tFieldFileName = tPrefix + "projects/HMR/test/data/hmr_Field_IO_Test.hdf5";

                std::shared_ptr< Field > tField = tHMR.load_field_from_file( "Circle", tFieldFileName,
                                                                             tLagrangeMeshIndex, tLBSplineMeshIndex );

                //            tHMR.save_to_exodus( tLagrangeMeshIndex, "HMR_Field_IO_Test_new.exo" );

                //------------------------------------------------------------------------------
                //   Test error
                //------------------------------------------------------------------------------

                auto tOutputMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                // calculate exact value
                auto tExact = tOutputMesh->create_field( "Exact", 0 );

                tExact->evaluate_scalar_function( LevelSetFunction );

                // determine coefficient of determination
                real tR2 = r2( tExact->get_node_values(), tField->get_node_values() );

                // perform test
                REQUIRE( tR2 > 0.99 );
            }
        }
    }
}
