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

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE( "Couette_Flow_Static",
        "[moris],[example],[fluid],[laminar]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = MORIS_EXAMPLE_INPUT_Couette_Flow_Static;

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    moris::exa::remove_test_outputs( { "Couette_Flow_Static.exo" } );
    REQUIRE( fn_WRK_Workflow_Main_Interface( argc, argv ) == 0 );
}
