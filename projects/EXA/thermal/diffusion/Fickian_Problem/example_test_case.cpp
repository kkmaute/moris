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

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------

// test case index
uint gTestCaseIndex;

// global variable for interpolation order
uint gInterpolationOrder;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------

extern "C" void
check_linear_results( moris::mtk::Exodus_IO_Helper &aExoIO, uint aNodeId )
{
    if ( gPrintReferenceValues )
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 2, 19 ) << '\n';

        // value of IQI at reference time step
        // std::cout << "IQI value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { { 3.920e-03 }, { 0.0 } };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );

    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.0e+02;

    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value() - tReferenceTime ) / tReferenceTime );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 3.442925667877556e+02;

    real tRelTempDifference = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 19 ) - tReferenceTemperature ) / tReferenceTemperature );

    // FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE( tRelTempDifference < 1.0e-5 );
}

//---------------------------------------------------------------

extern "C" void
check_quadratic_results( moris::mtk::Exodus_IO_Helper &aExoIO, uint aNodeId )
{
    if ( gPrintReferenceValues )
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 2, 19 ) << '\n';

        // value of IQI at reference time step
        // std::cout << "IQI value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { { 3.920e-03 }, { 0.0 } };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );

    // FIXME: fix global variable problem for quadratic
    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.0e+02;

    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value() - tReferenceTime ) / tReferenceTime );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 3.442923942621712e+02;

    real tRelTempDifference = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 19 ) - tReferenceTemperature ) / tReferenceTemperature );

    // FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE( tRelTempDifference < 1.0e-5 );
}

//---------------------------------------------------------------

extern "C" void
check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( "Fick_Problem_0.exo", 0, false, false );

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if ( gPrintReferenceValues )
    {
        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';
    }
    else
    {
        REQUIRE( tNumDims == 2 );
        REQUIRE( tNumNodes == 602 );
        REQUIRE( tNumElems == 300 );
    }

    // check results
    uint tNodeId = 14;

    check_linear_results( tExoIO, tNodeId );
}

//---------------------------------------------------------------

extern "C" void
check_quadratic_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( "Fick_Problem_1.exo", 0, false, false );

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if ( gPrintReferenceValues )
    {
        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';
    }
    else
    {
        REQUIRE( tNumDims == 2 );
        REQUIRE( tNumNodes == 1803 );
        REQUIRE( tNumElems == 300 );
    }

    // check results
    uint tNodeId = 39;

    check_quadratic_results( tExoIO, tNodeId );
}

//---------------------------------------------------------------

TEST_CASE( "Fickian_Problem_Linear",
        "[moris],[example],[thermal],[diffusion]" )
{
    // temporary
    gLogger.initialize( "Log.log" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Fick_Problem.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set test case index
    gTestCaseIndex = 0;

    // set interpolation order
    gInterpolationOrder = 1;

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            check_linear_results_serial();
            break;
        }
        default:
        {
            MORIS_ERROR( false, "This 2D Example can only be run in serial." );
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Fickian_Problem_Quadratic",
        "[moris],[example],[thermal],[diffusion]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Fick_Problem.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 2;

    // set test case index
    gTestCaseIndex = 1;

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            check_quadratic_results_serial();
            break;
        }
        default:
        {
            MORIS_ERROR( false, "This 2D Example can only be run in serial." );
        }
    }
}
