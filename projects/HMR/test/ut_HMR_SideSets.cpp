/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_HMR_SideSets.cpp
 *
 */

#include <catch.hpp>
#include "paths.hpp"

//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

//------------------------------------------------------------------------------
// from LINALG
#include "cl_Matrix.hpp"
#include "HDF5_Tools.hpp"
#include "op_equal_equal.hpp"
#include "fn_all_true.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"

//------------------------------------------------------------------------------

namespace moris::hmr
{

    TEST_CASE( "HMR_SideSets", "[moris],[mesh],[hmr],[hmr_side_set]" )
    {
        // get root from environment
        std::string tMorisRoot = get_base_moris_dir();

        //------------------------------------------------------------------------------
        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {

            SECTION( "HMR_SideSets 2D" )
            {
                // determine path for object file
                std::string tHdf5FilePath = tMorisRoot + "/projects/HMR/test/data/hmr_sideset_test_2d.hdf5";

                //------------------------------------------------------------------------------

                uint tLagrangeMeshIndex = 0;

                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 6 } } );
                tParameters.set_domain_dimensions( { { 4 }, { 6 } } );
                tParameters.set_domain_offset( { { 0.0 }, { 0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { 1 } );
                tParameters.set_lagrange_patterns( { 0 } );

                tParameters.set_bspline_orders( { 1 } );
                tParameters.set_bspline_patterns( { 0 } );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------

                HMR tHMR( tParameters );

                //------------------------------------------------------------------------------
                //    create refinement pattern
                //------------------------------------------------------------------------------

                std::shared_ptr< Database > tDatabase = tHMR.get_database();

                tDatabase->set_activation_pattern( tLagrangeMeshIndex );

                for ( uint tLevel = 0; tLevel < 4; ++tLevel )
                {
                    // flag first element
                    tDatabase->flag_element( 0 );

                    // flag last element
                    tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc() - 1 );

                    // manually refine, do not reset pattern
                    tDatabase->perform_refinement( 0, false );
                }

                // finish mesh
                tHMR.finalize();

                //------------------------------------------------------------------------------
                //    create Mesh and variables
                //------------------------------------------------------------------------------

                // create MTK mesh
                std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                Matrix< IndexMat > tElements;
                Matrix< IndexMat > tElementsSolution;
                Matrix< IndexMat > tOrds;

                //------------------------------------------------------------------------------
                //    write solution ( uncomment this if you want to recreate solution files )
                //------------------------------------------------------------------------------

                //                  // create file
                //                  hid_t tFileID = create_hdf5_file( tHdf5FilePath );
                //
                //                  // error handler
                //                  herr_t tStatus = 0;
                //
                //                  for( uint s=1; s<=4; ++s )
                //                  {
                //                      // create label
                //                      std::string tSetLabel = "SideSet_" + std::to_string( s );
                //
                //                      // ask mesh for sideset
                //                      tMesh->get_sideset_elems_loc_inds_and_ords(
                //                              tSetLabel, tElements, tOrds );
                //
                //                      print(tElements,"tElements");
                //
                //                      // save data
                //                      save_matrix_to_hdf5_file( tFileID, tSetLabel, tElements, tStatus );
                //                  }
                //
                //
                //                  // close file
                //                  close_hdf5_file( tFileID );
                //
                //                  // save exodus file for visual inspection
                //                  tHMR.save_to_exodus( tLagrangeMeshIndex, "Mesh.exo" );

                //------------------------------------------------------------------------------
                //    open solution
                //------------------------------------------------------------------------------

                // create file
                hid_t tFileID = open_hdf5_file( tHdf5FilePath );

                // error handler
                herr_t tStatus = 0;

                // loop over all sidesets
                for ( uint s = 1; s <= 4; ++s )
                {
                    // create label
                    std::string tSetLabel = "SideSet_" + std::to_string( s );

                    // ask mesh for sideset
                    tMesh->get_sideset_elems_loc_inds_and_ords( tSetLabel, tElements, tOrds );

                    // read solution from file
                    load_matrix_from_hdf5_file( tFileID, tSetLabel, tElementsSolution, tStatus );

                    // only test if solution is not empty ( if a sideset exists for this proc )
                    if ( tElementsSolution.length() > 0 )
                    {
                        // compare result
                        REQUIRE( all_true( tElements == tElementsSolution ) );
                    }
                }

                // close file
                close_hdf5_file( tFileID );
            }
                //------------------------------------------------------------------------------

            SECTION( "HMR_SideSets 3D" )
            {
                // determine path for object file
                std::string tHdf5FilePath = tMorisRoot + "/projects/HMR/test/data/hmr_sideset_test_3d.hdf5";

                //------------------------------------------------------------------------------

                uint tLagrangeMeshIndex = 0;

                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 4 }, { 6 }, { 10 } } );
                tParameters.set_domain_dimensions( { { 4 }, { 6 }, { 10 } } );
                tParameters.set_domain_offset( { { 0.0 }, { 0.0 }, { 0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { 1 } );
                tParameters.set_lagrange_patterns( { 0 } );

                tParameters.set_bspline_orders( { 1 } );
                tParameters.set_bspline_patterns( { 0 } );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------

                HMR tHMR( tParameters );

                //------------------------------------------------------------------------------
                //    create refinement pattern
                //------------------------------------------------------------------------------

                std::shared_ptr< Database > tDatabase = tHMR.get_database();

                tDatabase->set_activation_pattern( tLagrangeMeshIndex );

                for ( uint tLevel = 0; tLevel < 4; ++tLevel )
                {
                    // flag first element
                    tDatabase->flag_element( 0 );

                    // flag last element
                    tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc() - 1 );

                    // manually refine, do not reset pattern
                    tDatabase->perform_refinement( tLagrangeMeshIndex, false );
                }

                // finish mesh
                tHMR.finalize();

                //------------------------------------------------------------------------------
                //    create Mesh and variables
                //------------------------------------------------------------------------------

                // create MTK mesh
                std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                Matrix< IndexMat > tElements;
                Matrix< IndexMat > tElementsSolution;
                Matrix< IndexMat > tOrds;

                //------------------------------------------------------------------------------
                //    write solution ( uncomment this if you want to recreate solution files )
                //------------------------------------------------------------------------------

                //          // create file
                //          hid_t tFileID = create_hdf5_file( tHdf5FilePath );
                //          //error handler
                //          herr_t tStatus = 0;
                //
                //          for( uint s=1; s<=6; ++s )
                //          {
                //              // create label
                //              std::string tSetLabel = "SideSet_" + std::to_string( s );
                //
                //              // ask mesh for sideset
                //              tMesh->get_sideset_elems_loc_inds_and_ords(
                //                      tSetLabel, tElements, tOrds );
                //
                //              // save data
                //              save_matrix_to_hdf5_file( tFileID, tSetLabel, tElements, tStatus );
                //          }
                //
                //          // close file
                //          close_hdf5_file( tFileID );

                //------------------------------------------------------------------------------
                //    open solution
                //------------------------------------------------------------------------------

                // create file
                hid_t tFileID = open_hdf5_file( tHdf5FilePath );

                // error handler
                herr_t tStatus = 0;

                // loop over all sidesets
                for ( uint s = 1; s <= 6; ++s )
                {
                    // create label
                    std::string tSetLabel = "SideSet_" + std::to_string( s );

                    // ask mesh for sideset
                    tMesh->get_sideset_elems_loc_inds_and_ords( tSetLabel, tElements, tOrds );

                    // read solution from file
                    load_matrix_from_hdf5_file( tFileID, tSetLabel, tElementsSolution, tStatus );

                    // only test if solution is not empty ( if a sideset exists for this proc )
                    if ( tElementsSolution.length() > 0 )
                    {
                        // compare result
                        //                  REQUIRE( all_true( tElements == tElementsSolution ) );
                    }
                }

                // close file
                close_hdf5_file( tFileID );
            }
        }

        //------------------------------------------------------------------------------
    }

    TEST_CASE( "HMR_SideSets_numbered_aura", "[moris],[mesh],[hmr],[hmr_side_set_numbered_aura]" )
    {
        // get root from environment
        std::string tMorisRoot = get_base_moris_dir();

        //------------------------------------------------------------------------------
        if ( par_size() == 1 || par_size() == 2 || par_size() == 4 )
        {

            SECTION( "HMR_SideSets 2D" )
            {
                // determine path for object file
                std::string tHdf5FilePath = tMorisRoot + "/projects/HMR/test/data/hmr_sideset_test_numbered_aura_2d.hdf5";

                //------------------------------------------------------------------------------

                uint tLagrangeMeshIndex = 0;

                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 6 }, { 6 } } );
                tParameters.set_domain_dimensions( { { 4 }, { 6 } } );
                tParameters.set_domain_offset( { { 0.0 }, { 0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 } } );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { 1 } );
                tParameters.set_lagrange_patterns( { 0 } );

                tParameters.set_bspline_orders( { 1 } );
                tParameters.set_bspline_patterns( { 0 } );

                tParameters.set_number_aura( true );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------

                HMR tHMR( tParameters );

                //------------------------------------------------------------------------------
                //    create refinement pattern
                //------------------------------------------------------------------------------

                std::shared_ptr< Database > tDatabase = tHMR.get_database();

                tDatabase->set_activation_pattern( tLagrangeMeshIndex );

                for ( uint tLevel = 0; tLevel < 4; ++tLevel )
                {
                    // flag first element
                    tDatabase->flag_element( 0 );

                    // flag last element
                    tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc() - 1 );

                    // manually refine, do not reset pattern
                    tDatabase->perform_refinement( 0, false );
                }

                // finish mesh
                tHMR.finalize();

                //------------------------------------------------------------------------------
                //    create Mesh and variables
                //------------------------------------------------------------------------------

                // create MTK mesh
                std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                Matrix< IndexMat > tElements;
                Matrix< IndexMat > tElementsSolution;
                Matrix< IndexMat > tOrds;

                //------------------------------------------------------------------------------
                //    write solution ( uncomment this if you want to recreate solution files )
                //------------------------------------------------------------------------------

                //                  // create file
                //                  hid_t tFileID = create_hdf5_file( tHdf5FilePath );
                //
                //                  // error handler
                //                  herr_t tStatus = 0;
                //
                //                  for( uint s=1; s<=4; ++s )
                //                  {
                //                      // create label
                //                      std::string tSetLabel = "SideSet_" + std::to_string( s );
                //
                //                      // ask mesh for sideset
                //                      tMesh->get_sideset_elems_loc_inds_and_ords(
                //                              tSetLabel, tElements, tOrds );
                //
                //                      // save data
                //                      save_matrix_to_hdf5_file( tFileID, tSetLabel, tElements, tStatus );
                //                  }
                //
                //
                //                  // close file
                //                  close_hdf5_file( tFileID );
                //
                //                  // save exodus file for visual inspection
                //                  tHMR.save_to_exodus( tLagrangeMeshIndex, "Mesh.exo" );

                //------------------------------------------------------------------------------
                //    open solution
                //------------------------------------------------------------------------------
                //
                // create file
                hid_t tFileID = open_hdf5_file( tHdf5FilePath );

                // error handler
                herr_t tStatus = 0;

                // loop over all sidesets
                for ( uint s = 1; s <= 4; ++s )
                {
                    // create label
                    std::string tSetLabel = "SideSet_" + std::to_string( s );

                    // ask mesh for sideset
                    tMesh->get_sideset_elems_loc_inds_and_ords( tSetLabel, tElements, tOrds );

                    // read solution from file
                    load_matrix_from_hdf5_file( tFileID, tSetLabel, tElementsSolution, tStatus );

                    // only test if solution is not empty ( if a sideset exists for this proc )
                    if ( tElementsSolution.length() > 0 )
                    {
                        // compare result
                        REQUIRE( all_true( tElements == tElementsSolution ) );
                    }
                }

                // close file
                close_hdf5_file( tFileID );
            }
                //------------------------------------------------------------------------------

            SECTION( "HMR_SideSets 3D" )
            {
                // determine path for object file
                std::string tHdf5FilePath = tMorisRoot + "/projects/HMR/test/data/hmr_sideset_test_numbered_aura_3d.hdf5";

                //------------------------------------------------------------------------------

                uint tLagrangeMeshIndex = 0;

                Parameters tParameters;

                tParameters.set_number_of_elements_per_dimension( { { 6 }, { 6 }, { 10 } } );
                tParameters.set_domain_dimensions( { { 4 }, { 6 }, { 10 } } );
                tParameters.set_domain_offset( { { 0.0 }, { 0.0 }, { 0.0 } } );
                tParameters.set_bspline_truncation( true );

                tParameters.set_side_sets( { { 1 }, { 2 }, { 3 }, { 4 }, { 5 }, { 6 } } );

                tParameters.set_output_meshes( { { { 0 } } } );

                tParameters.set_lagrange_orders( { 1 } );
                tParameters.set_lagrange_patterns( { 0 } );

                tParameters.set_bspline_orders( { 1 } );
                tParameters.set_bspline_patterns( { 0 } );

                Vector< Matrix< DDSMat > > tLagrangeToBSplineMesh( 1 );
                tLagrangeToBSplineMesh( 0 ) = { { 0 } };

                tParameters.set_number_aura( true );

                tParameters.set_lagrange_to_bspline_mesh( tLagrangeToBSplineMesh );

                //------------------------------------------------------------------------------

                HMR tHMR( tParameters );

                //------------------------------------------------------------------------------
                //    create refinement pattern
                //------------------------------------------------------------------------------

                std::shared_ptr< Database > tDatabase = tHMR.get_database();

                tDatabase->set_activation_pattern( tLagrangeMeshIndex );

                for ( uint tLevel = 0; tLevel < 4; ++tLevel )
                {
                    // flag first element
                    tDatabase->flag_element( 0 );

                    // flag last element
                    tDatabase->flag_element( tDatabase->get_number_of_elements_on_proc() - 1 );

                    // manually refine, do not reset pattern
                    tDatabase->perform_refinement( tLagrangeMeshIndex, false );
                }

                // finish mesh
                tHMR.finalize();

                //------------------------------------------------------------------------------
                //    create Mesh and variables
                //------------------------------------------------------------------------------

                // create MTK mesh
                std::shared_ptr< Mesh > tMesh = tHMR.create_mesh( tLagrangeMeshIndex );

                Matrix< IndexMat > tElements;
                Matrix< IndexMat > tElementsSolution;
                Matrix< IndexMat > tOrds;

                //------------------------------------------------------------------------------
                //    write solution ( uncomment this if you want to recreate solution files )
                //------------------------------------------------------------------------------

                //          // create file
                //          hid_t tFileID = create_hdf5_file( tHdf5FilePath );
                //          //error handler
                //          herr_t tStatus = 0;
                //
                //          for( uint s=1; s<=6; ++s )
                //          {
                //              // create label
                //              std::string tSetLabel = "SideSet_" + std::to_string( s );
                //
                //              // ask mesh for sideset
                //              tMesh->get_sideset_elems_loc_inds_and_ords(
                //                      tSetLabel, tElements, tOrds );
                //
                //              // save data
                //              save_matrix_to_hdf5_file( tFileID, tSetLabel, tElements, tStatus );
                //          }
                //
                //          // close file
                //          close_hdf5_file( tFileID );

                //------------------------------------------------------------------------------
                //    open solution
                //------------------------------------------------------------------------------

                // create file
                hid_t tFileID = open_hdf5_file( tHdf5FilePath );

                // error handler
                herr_t tStatus = 0;

                // loop over all sidesets
                for ( uint s = 1; s <= 6; ++s )
                {
                    // create label
                    std::string tSetLabel = "SideSet_" + std::to_string( s );

                    // ask mesh for sideset
                    tMesh->get_sideset_elems_loc_inds_and_ords( tSetLabel, tElements, tOrds );

                    // read solution from file
                    load_matrix_from_hdf5_file( tFileID, tSetLabel, tElementsSolution, tStatus );

                    // only test if solution is not empty ( if a sideset exists for this proc )
                    if ( tElementsSolution.length() > 0 )
                    {
                        // compare result
                        //                  REQUIRE( all_true( tElements == tElementsSolution ) );
                    }
                }

                // close file
                close_hdf5_file( tFileID );
            }
        }

        //------------------------------------------------------------------------------
    }
}
