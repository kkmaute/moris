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
        const std::string& aExoFileName,
        const std::string& aLabel )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results for stage %s on %i processor.", aLabel.c_str(), par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    // reference node IDs to check
    Vector< uint > tNodeIDs = { 0, 100, 200, 300, 400, 500, 600, 700, 800 };

    std::cout << "Checking: " << aLabel << '\n';

    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if ( gPrintReferenceValues )
    {

        std::cout << "uint tReferenceNumDims     = " << tNumDims << ";\n";
        std::cout << "uint tReferenceNumNodes    = " << tNumNodes << ";\n";
        std::cout << "uint tReferenceNumElements = " << tNumElems << ";\n";

        std::cout << "Vector< real > tLevelSetValues = {\n";

        for ( uint iCheckNode = 0; iCheckNode < tNodeIDs.size(); iCheckNode++ )
        {
            if ( iCheckNode > 0 )
            {
                std::cout << ",\n";
            }
            // get nodal level set value
            std::cout << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tNodeIDs( iCheckNode ), 2, 0 );
        }
        std::cout << "};\n";

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    uint tReferenceNumDims     = 2;
    uint tReferenceNumNodes    = 12103;
    uint tReferenceNumElements = 9997;

    // define reference nodal values
    Vector< real > tLevelSetValues = {
        -1.103832737029360e-04,
        7.297122490942622e+00,
        1.158782219523812e+01,
        1.333500196238856e+01,
        1.679078682912798e-02,
        -1.965913331381824e-01,
        -6.847209299360845e+00,
        -1.052316058716752e+01,
        -1.329570299466918e+01
    };

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims,
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims ) / tReferenceNumDims * 100.0 ) );

    MORIS_LOG_INFO( "Check number of nodes: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes,
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes ) / tReferenceNumNodes * 100.0 ) );

    MORIS_LOG_INFO( "Check number of elements: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElements,
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElements ) / tReferenceNumElements * 100.0 ) );

    // Check field values
    for ( uint iCheckNode = 0; iCheckNode < 9; iCheckNode++ )
    {
        real tValue = tExoIO.get_nodal_field_value( tNodeIDs( iCheckNode ), 2, 0 );

        MORIS_LOG_INFO( "Check level set value at node %d: reference %12.5e, actual %12.5e, percent  error %12.5e.",
                tNodeIDs( iCheckNode ),
                tLevelSetValues( iCheckNode ),
                tValue,
                std::abs( ( tValue - tLevelSetValues( iCheckNode ) ) / tLevelSetValues( iCheckNode ) * 100.0 ) );

        CHECK( tValue == Approx( tLevelSetValues( iCheckNode ) ).epsilon( 0.0001 ) );
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

    // perform check on initial file
    check_results( "Levelset_Boxbeam_Create_File.exo.e-s.0011", "Create" );
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

    // perform check on restart file
    check_results( "Levelset_Boxbeam_Restart.exo.e-s.0011", "Restart" );
}
