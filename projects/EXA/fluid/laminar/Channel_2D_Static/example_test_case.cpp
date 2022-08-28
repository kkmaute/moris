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

#include "cl_Logger.hpp" // MRS/IOS/src

//---------------------------------------------------------------

// global variable for inlet BC
bool gInletVelocityBCFlag = true;
bool gInletPressureBCFlag = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("Channel_2D_Static_Inlet_Velocity",
        "[moris],[example],[fluid],[laminar]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_2D_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    //REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Channel_2D_Static_Inlet_Pressure",
        "[moris],[example],[fluid],[laminar]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_2D_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    //REQUIRE( tRet ==  0 );
}

