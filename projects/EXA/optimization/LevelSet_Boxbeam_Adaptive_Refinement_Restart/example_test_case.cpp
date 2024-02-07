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

#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"

using namespace moris;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        std::string aExoFileName,
        uint        aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        return;
    }

    // Reference nodal values to check
    Vector< uint > tNodeIDs = { 0, 100, 200, 300, 400, 500, 600, 700, 800 };
    Vector< real > tLevelSetValues = { -0.000135399, 7.40848, 14.3591, 15.4603, 0.861916, -0.499474, -7.59165, -13.1653, -15.2731 };

    // Check field values
    for ( uint iCheckNode = 0; iCheckNode < 9; iCheckNode++ )
    {
        CHECK( tExoIO.get_nodal_field_value( tNodeIDs( iCheckNode ), 2, 0 ) == Approx( tLevelSetValues( iCheckNode ) ).epsilon( 0.0001 ) );
    }
}

//---------------------------------------------------------------

TEST_CASE( "Leveset Boxbeam Create File",
        "[moris],[example],[optimization],[levelset_boxbeam_create_file]" )
{
    // remove files from previous test runs
    // FIXME: should be made independent of OS; note std::remove does not take wild cards
    if ( par_rank() == 0 )
    {
        MORIS_ERROR( std::system( "rm -f *exo*" ) == 0, "Leveset Boxbeam Create File - removing *exo* files failed" );
        MORIS_ERROR( std::system( "rm -f *hdf5*" ) == 0, "Leveset Boxbeam Create File - removing *hdf5* files failed" );
    }

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Create_File.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );
}

//---------------------------------------------------------------

TEST_CASE( "Leveset Boxbeam Restart",
        "[moris],[example],[optimization],[levelset_boxbeam_restart]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Restart.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set test case index
    uint tTestCaseIndex = 0;

    // perform check for Test Case 0
    check_results( "Levelset_Boxbeam_Restart.exo.e-s.0016", tTestCaseIndex );
}
