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

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

// initialize various solver parameters as global variables

bool gHaveStaggeredFA = false;
bool gHaveStaggeredSA = false;  // only relevant if Fwd. analysis also staggered
bool gUseMixedTimeElements = false;
bool gUseBelosWithILUT = false;

//---------------------------------------------------------------

TEST_CASE("Standard_Monolithic",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Standard_Monolithic]")
{
    // change parameters
    gHaveStaggeredFA = false;
    gHaveStaggeredSA = false;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

// FIXME: This solver configuration returns an error in the matrix assembly due to the mixed time elements
// TEST_CASE("Monolithic_Mixed_Time_Elements",
//         "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Monolithic_Mixed_Time_Elements]")
// {
//     // change parameters
//     gHaveStaggeredFA = false;
//     gHaveStaggeredSA = false;
//     gUseMixedTimeElements = true;
//     gUseBelosWithILUT = false;
//
//     // define command line call
//     int argc = 2;
//
//     char tString1[] = "";
//     char tString2[] = "Solver_Examples.so";
//
//     char * argv[2] = {tString1,tString2};
//
//     // call to performance manager main interface
//     int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
//
//     // catch test statements should follow
//     REQUIRE( tRet ==  0 );
// }

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_Monolithic_SA",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_Monolithic_SA]")
{
    // change parameters
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = false;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_and_SA",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_and_SA]")
{
    // change parameters
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = true;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_and_SA_Mixed_Time_Elements",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_and_SA_Mixed_Time_Elements]")
{
    // change parameters
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = true;
    gUseMixedTimeElements = true;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

