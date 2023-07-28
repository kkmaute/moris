/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 */

#include <catch.hpp>

#include "cl_Logger.hpp"    // MRS/IOS/src
#include "HDF5_Tools.hpp"

using namespace moris;

// global variable to define test cases
uint tGeoModel;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results( uint aTestCaseIndex, std::string aHDF5FileName)
{
    MORIS_LOG_INFO( "" );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processors.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( "" );

    // Declare sensitivity matrices for comparison
    Matrix< DDRMat > tObjectiveAnalytical;
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tObjectiveFD;
    Matrix< DDRMat > tConstraintsFD;

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( aHDF5FileName );
    herr_t tStatus = 0;

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytical, tStatus );
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );
    REQUIRE( tObjectiveAnalytical.length() == tConstraintsAnalytical.length() );    // one objective and one constraint for this problem only

    // Read FD sensitivities and compare
    Cell< std::string > tFDTypes = { "fd_forward", "fd_backward", "fd_central" };
    for ( uint tFDIndex = 0; tFDIndex < tFDTypes.size(); tFDIndex++ )
    {
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tObjectiveFD, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tConstraintsFD, tStatus );

        REQUIRE( tObjectiveAnalytical.length() == tObjectiveFD.length() );
        REQUIRE( tConstraintsAnalytical.length() == tConstraintsFD.length() );

        for ( uint tADVIndex = 0; tADVIndex < tObjectiveAnalytical.length(); tADVIndex++ )
        {
            MORIS_LOG_INFO( "Check derivative of objective  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tObjectiveAnalytical( tADVIndex ),
                    tFDTypes( tFDIndex ).c_str(),
                    tObjectiveFD( tADVIndex ),
                    100 * std::abs( ( tObjectiveAnalytical( tADVIndex ) - tObjectiveFD( tADVIndex ) ) / ( tObjectiveFD( tADVIndex ) + MORIS_REAL_EPS ) ) );

            MORIS_LOG_INFO( "Check derivative of constraint wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tConstraintsAnalytical( tADVIndex ),
                    tFDTypes( tFDIndex ).c_str(),
                    tConstraintsFD( tADVIndex ),
                    100 * std::abs( ( tConstraintsAnalytical( tADVIndex ) - tConstraintsFD( tADVIndex ) ) / ( tConstraintsFD( tADVIndex ) + MORIS_REAL_EPS ) ) );

            CHECK( tObjectiveAnalytical( tADVIndex ) == Approx( tObjectiveFD( tADVIndex ) ).margin( MORIS_REAL_EPS ) );
            CHECK( tConstraintsAnalytical( tADVIndex ) == Approx( tConstraintsFD( tADVIndex ) ).margin( MORIS_REAL_EPS ) );
        }
    }

    // close file
    close_hdf5_file( tFileID );
}

// ---------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Shape_Sensitivity_Bspline_2D",
        "[moris],[example],[optimization],[sweep]" )
{
    // remove files from previous test runs
    // FIXME: should be made independent of OS; note std::remove does not take wild cards
    if ( par_rank() == 0 )
    {
        MORIS_ERROR( std::system( "rm -f *exo*" ) == 0, "Shape_Sensitivity_Bspline - removing *exo* files failed" );
        MORIS_ERROR( std::system( "rm -f *hdf5*" ) == 0, "Shape_Sensitivity_Bspline - removing *hdf5* files failed" );
    }

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Shape_Sensitivity_Bspline_2D.so";

    char* argv[ 2 ] = { tString1, tString2 };

    if ( par_size() == 1 )
    {
        // loop over all test configurations
        for ( tGeoModel = 0; tGeoModel < 7; ++tGeoModel )
        {
            // call to performance manager main interface
            int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

            // catch test statements should follow
            REQUIRE( tRet == 0 );

            // check results
            check_results( tGeoModel, "shape_opt_test_2D.hdf5" );
        }
    }
}

TEST_CASE( "Shape_Sensitivity_Bspline_3D",
        "[moris],[example],[optimization],[sweep]" )
{
    // remove files from previous test runs
    // FIXME: should be made independent of OS; note std::remove does not take wild cards
    if ( par_rank() == 0 )
    {
        MORIS_ERROR( std::system( "rm -f *exo*" ) == 0, "Shape_Sensitivity_Bspline - removing *exo* files failed" );
        MORIS_ERROR( std::system( "rm -f *hdf5*" ) == 0, "Shape_Sensitivity_Bspline - removing *hdf5* files failed" );
    }

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Shape_Sensitivity_Bspline_3D.so";

    char* argv[ 2 ] = { tString1, tString2 };

    if ( par_size() == 1 )
    {
        // loop over all test configurations
        for ( tGeoModel = 0; tGeoModel < 7; ++tGeoModel )
        {
            // call to performance manager main interface
            int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

            // catch test statements should follow
            REQUIRE( tRet == 0 );

            // check results
            check_results( tGeoModel, "shape_opt_test_3D.hdf5" );
        }
    }
}
