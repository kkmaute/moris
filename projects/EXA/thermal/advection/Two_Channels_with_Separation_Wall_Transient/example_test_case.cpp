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

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("Two_Channels_with_Separation_Wall_Transient",
        "[moris],[example],[thermal],[advection]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./bin/Two_Channels_with_Separation_Wall_Transient.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    //REQUIRE( tRet ==  0 );

//    // open and query exodus output file (set verbose to true to get basic mesh information)
//    moris::mtk::Exodus_IO_Helper tExoIO("Two_Channels_with_Separation_Wall_Transient.exo",0,false);
//
//    // check dimension, number of nodes and number of elements
//    uint tNumDims = tExoIO.get_number_of_dimensions();
//
//    REQUIRE( tNumDims                          ==  2    );
//    REQUIRE( tExoIO.get_number_of_nodes()      ==  1081 );
//    REQUIRE( tExoIO.get_number_of_elements()   ==  1100 );
//
//    // define reference coordinates for node 488
//    // moris::print( tExoIO.get_nodal_coordinate( 488 ), "");
//    Matrix< DDRMat > tReferenceCoordinate = { {+8.509090909090911e-01},{2.309090909090910e-01} };
//
//    // check nodal coordinates
//    real tRelDiffNorm = moris::norm( tExoIO.get_nodal_coordinate( 488 ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);
//
//    REQUIRE( tRelDiffNorm <  1.0e-8 );
//
//    // check time value for time step index 0
//    // std::cout << std::scientific << std::setprecision(15) << tExoIO.get_time_value() << std::endl;
//    real tReferenceTime = 1.000000000000000e+00;
//
//    real tRelTimeDifference = std::abs( ( tExoIO.get_time_value( ) - tReferenceTime) / tReferenceTime );
//
//    REQUIRE( tRelTimeDifference <  1.0e-8 );
//
//    // check temperature at node 488 in first time step (temperature is 3rd nodal field, first time step has index 0)
//    // std::cout << std::scientific << std::setprecision(15) << tExoIO.get_nodal_field_value( 488, 2, 0 ) << std::endl;
//    real tReferenceTemperature = 9.343574182982935e+04;
//
//    real tRelTempDifference = std::abs( ( tExoIO.get_nodal_field_value( 488, 2, 0 ) - tReferenceTemperature ) / tReferenceTemperature );
//    REQUIRE(  tRelTempDifference < 1.0e-8);
//
//    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
//    // std::cout << std::scientific << std::setprecision(15) << tExoIO.get_global_variable(0, 0 ) << std::endl;
//    real tReferenceIQI = 8.068274679000535e+04;
//
//    real tRelIQIDifference = std::abs( ( tExoIO.get_global_variable(0, 0 ) - tReferenceIQI ) / tReferenceIQI );
//    REQUIRE(  tRelIQIDifference < 1.0e-8);
}

