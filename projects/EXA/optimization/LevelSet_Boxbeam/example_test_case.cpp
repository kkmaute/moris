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

    if ( true )
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

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2 };
    Vector< uint > tReferenceNumNodes = { 2451 };
    Vector< uint > tReferenceNumElems = { 1875 };

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

TEST_CASE( "Leveset Boxbeam",
        "[moris],[example],[optimization],[levelset_boxbeam]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Levelset_Boxbeam.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // Read Exodus file
    // moris::mtk::Exodus_IO_Helper tExoIO("levelset_boxbeam.exo", 0, false, false);

    // Checks
}

// TEST_CASE("LevelSet_BoxBeam_LinQuad",
//         "[moris],[example],[optimization],[levelset_boxbeam]")
//{
//     if( moris::par_size() == 1 )
//     {
//         // define command line call
//         int argc = 2;
//
//         char tString1[] = "";
//         char tString2[] = "Levelset_Boxbeam_Lin_Quad_LevelSet.so";
//
//         char * argv[2] = {tString1,tString2};
//
//         // call to performance manager main interface
//         int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
//
//         // catch test statements should follow
//         REQUIRE( tRet ==  0 );
//
//         // set test case index
//         uint tTestCaseIndex = 0;
//
//         // Read Exodus file
//         //moris::mtk::Exodus_IO_Helper tExoIO("levelset_boxbeam.exo", 0, false, false);
//
//         // Checks
//         // perform check for Test Case 0
//         check_results("Levelset_Boxbeam_Lin_Quad_LevelSet.exo.e-s.0002",tTestCaseIndex);
//     }
// }
