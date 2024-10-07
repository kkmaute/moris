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
check_results_serial(
        const std::string& aExoFileName,
        uint               aTestCaseIndex )
{

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
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
    Vector< uint > tReferenceNumNodes = { 4495, 17545 };
    Vector< uint > tReferenceNumElems = { 4182, 4236 };

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
        const std::string& aExoFileName,
        uint               aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << " on proc: " << par_rank() << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Rank " << par_rank() << " Number of dimensions: " << tNumDims << '\n';
        std::cout << "Rank " << par_rank() << " Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Rank " << par_rank() << " Number of elements  : " << tNumElems << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims = { 2, 2 };
    Vector< uint > tReferenceNumNodes;
    Vector< uint > tReferenceNumElems;

    if ( par_rank() == 0 )
    {
        tReferenceNumNodes = { 1276, 5553 };
        tReferenceNumElems = { 1157, 1320 };
    }
    if ( par_rank() == 1 )
    {
        tReferenceNumNodes = { 1276, 5553 };
        tReferenceNumElems = { 1157, 1320 };
    }
    if ( par_rank() == 2 )
    {
        tReferenceNumNodes = { 1267, 5553 };
        tReferenceNumElems = { 1148, 1320 };
    }
    if ( par_rank() == 3 )
    {
        tReferenceNumNodes = { 1267, 5553 };
        tReferenceNumElems = { 1148, 1320 };
    }

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    // std::cout<<"nodes: "<<tNumNodes<<" ele: "<<tNumElems<<" proc: "<<par_rank()<<std::endl;

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

TEST_CASE( "Leveset Boxbeam Adaptive Refinement",
        "[moris],[example],[optimization],[levelset_boxbeam_adaptive_refinement]" )
{
    // remove files from previous test runs
    // FIXME: should be made independent of OS; note std::remove does not take wild cards
    if ( par_rank() == 0 )
    {
        MORIS_ERROR( std::system( "rm -f *exo*" ) == 0, "Leveset Boxbeam Adaptive Refinement - removing *exo* files failed" );
        MORIS_ERROR( std::system( "rm -f *hdf5*" ) == 0, "Leveset Boxbeam Adaptive Refinement - removing *hdf5* files failed" );
    }

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Adaptive_Refinement.so";

    char* argv[ 2 ] = { tString1, tString2 };

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
            check_results_serial( "HMRLagrangeMesh.exo", tTestCaseIndex );
            break;
        }
        case 4:
        {
            // perform check
            check_results_parallel( "HMRLagrangeMesh.exo", tTestCaseIndex );

            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Leveset Boxbeam Adaptive Refinement Quadratic",
        "[moris],[example],[optimization],[levelset_boxbeam_adaptive_refinement_quadratic]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam_Adaptive_Refinement_Quadratic.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set test case index
    uint tTestCaseIndex = 1;

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            // perform check for Test Case 0
            check_results_serial( "HMRLagrangeMesh_Quadratic.exo", tTestCaseIndex );
            break;
        }
        case 4:
        {

            // perform check
            check_results_parallel( "HMRLagrangeMesh_Quadratic.exo", tTestCaseIndex );
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}
