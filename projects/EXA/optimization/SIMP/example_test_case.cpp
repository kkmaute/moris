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

TEST_CASE("SIMP",
        "[moris],[example],[optimization],[sweep]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "SIMP.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );

    // Read Exodus file
    moris::mtk::Exodus_IO_Helper tExoIO("SIMP.exo", 0, false, false);

    // Checks
    CHECK(tExoIO.get_nodal_field_value(   0, 2, 0 ) == Approx(0.419977));
    CHECK(tExoIO.get_nodal_field_value( 100, 2, 0 ) == Approx(0.419974));
    CHECK(tExoIO.get_nodal_field_value( 200, 2, 0 ) == Approx(0.380038));
    CHECK(tExoIO.get_nodal_field_value( 300, 2, 0 ) == Approx(0.397622));
    CHECK(tExoIO.get_nodal_field_value( 400, 2, 0 ) == Approx(0.380049));
    CHECK(tExoIO.get_nodal_field_value( 500, 2, 0 ) == Approx(0.380093));
    CHECK(tExoIO.get_nodal_field_value( 600, 2, 0 ) == Approx(0.380112));
    CHECK(tExoIO.get_nodal_field_value( 700, 2, 0 ) == Approx(0.419978));
    CHECK(tExoIO.get_nodal_field_value( 800, 2, 0 ) == Approx(0.419973));

}

