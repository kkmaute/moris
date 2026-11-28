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
#include "fn_EXA_Test_Utilities.hpp"

#include "cl_Logger.hpp"    // MRS/IOS/src

//---------------------------------------------------------------

// global variable for inlet BC
bool gInletVelocityBCFlag = true;
bool gInletPressureBCFlag = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE( "Channel_2D_Static_Inlet_Velocity",
        "[moris],[example],[fluid],[laminar]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = MORIS_EXAMPLE_INPUT_Channel_2D_Static;

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    moris::exa::remove_test_outputs( { "Channel_2D_Static.exo" } );
    REQUIRE( fn_WRK_Workflow_Main_Interface( argc, argv ) == 0 );
}

//---------------------------------------------------------------

TEST_CASE( "Channel_2D_Static_Inlet_Pressure",
        "[moris],[example],[fluid],[laminar]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = MORIS_EXAMPLE_INPUT_Channel_2D_Static;

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    moris::exa::remove_test_outputs( { "Channel_2D_Static.exo" } );
    REQUIRE( fn_WRK_Workflow_Main_Interface( argc, argv ) == 0 );
}
