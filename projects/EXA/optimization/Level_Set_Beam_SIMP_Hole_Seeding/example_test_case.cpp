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

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//------------------------------------------------------------------------------------

void test_pause()
{
    // communicate process ids
    pid_t tPId = getpid();

    moris::Matrix< moris::DDSMat > tPIdVec;
    allgather_scalar( tPId, tPIdVec );

    if ( moris::par_rank() == 0 )
    {
        // print process Ids
        for ( int i = 0; i < moris::par_size(); ++i )
        {
            fprintf( stderr, "Process Rank %d ID: %d\n", i, tPIdVec( i ) );
        }

        std::string dummy;

        std::cout << "Press enter to continue . . .\n"
                  << std::flush;
        std::getline( std::cin, dummy );
    }
}

//---------------------------------------------------------------

extern "C" void
check_results_serial(
        const std::string &aExoFileName,
        uint               aTestCaseIndex )
{

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( true )
    {
        std::cout << "Test case index: " << aTestCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2 };
    Vector< uint > tReferenceNumNodes = { 10477, 13815 };
    Vector< uint > tReferenceNumElems = { 8178, 7166 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );
}

extern "C" void
check_results_parallel(
        const std::string &aExoFileName,
        uint               aTestCaseIndex )
{

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( true )
    {
        std::cout << "Test case index: " << aTestCaseIndex << " on proc: " << par_rank() << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims = { 2, 2 };
    Vector< uint > tReferenceNumNodes;
    Vector< uint > tReferenceNumElems;

    if ( par_rank() == 0 )
    {
        tReferenceNumNodes = { 2739, 3488 };
        tReferenceNumElems = { 2039, 1768 };
    }
    if ( par_rank() == 1 )
    {
        tReferenceNumNodes = { 2738, 3492 };
        tReferenceNumElems = { 2041, 1771 };
    }
    if ( par_rank() == 2 )
    {
        tReferenceNumNodes = { 2560, 3379 };
        tReferenceNumElems = { 2050, 1746 };
    }
    if ( par_rank() == 3 )
    {
        tReferenceNumNodes = { 2574, 3630 };
        tReferenceNumElems = { 2047, 1884 };
    }

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );
}

//---------------------------------------------------------------

TEST_CASE( "Level_Set_Beam_SIMP_Hole_Seeding",
        "[moris],[example],[optimization],[Level_Set_Beam_SIMP_Hole_Seeding]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Level_Set_Beam_SIMP_Hole_Seeding.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // for debugging in parallel
    test_pause();

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set test case index
    uint tTestCaseIndex = 0;

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            // perform check for Test Case 0
            check_results_serial( "Level_Set_Beam_SIMP_Hole_Seeding.exo.e-s.0018", tTestCaseIndex );
            break;
        }
        case 4:
        {
            // perform check
            check_results_parallel( "Level_Set_Beam_SIMP_Hole_Seeding.exo.e-s.0018", tTestCaseIndex );

            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

//---------------------------------------------------------------
