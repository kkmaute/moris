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
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 2, 2 ) << '\n';

        // solution of reference point at reference time step
        std::cout << "x-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 3, 2 ) << '\n';

        // solution of reference point at reference time step
        std::cout << "y-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 4, 2 ) << '\n';

        // value of Volume-IQI at reference time step
        std::cout << "Volume-IQI value: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_global_variable( 0, 0 ) << '\n';

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { { 0.0 }, { 5.0e-02 } };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );
    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime     = 0.10;
    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value() - tReferenceTime ) / tReferenceTime );
    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 1.334186887938226e+00;
    real tRelTempDifference    = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 2 ) - tReferenceTemperature ) / tReferenceTemperature );
    REQUIRE( tRelTempDifference < 1.0e-8 );

    // check x-displacement at node aNodeId in first time step
    real tReferenceXDisplacement = -8.310495618961509e-04;
    real tRelXDispDifference     = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 3, 2 ) - tReferenceXDisplacement ) / tReferenceXDisplacement );
    REQUIRE( tRelXDispDifference < 1.0e-6 );

    // check x-displacement at node aNodeId in first time step
    real tReferenceYDisplacement = 5.243129489358523e-03;
    real tRelYDispDifference     = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 4, 2 ) - tReferenceYDisplacement ) / tReferenceYDisplacement );
    REQUIRE( tRelYDispDifference < 1.0e-6 );

    // check volume-IQI at node aNodeId in first time step
    real tReferenceVolume     = 5.00e-05;
    real tRelVolumeDifference = std::abs( ( aExoIO.get_global_variable( 0, 0 ) - tReferenceVolume ) / tReferenceVolume );
    REQUIRE( tRelVolumeDifference < 1.0e-8 );
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
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 2, 2 ) << '\n';

        // solution of reference point at reference time step
        std::cout << "x-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 3, 2 ) << '\n';

        // solution of reference point at reference time step
        std::cout << "y-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 4, 2 ) << '\n';

        // value of Volume-IQI at reference time step
        std::cout << "Volume-IQI value: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_global_variable( 0, 0 ) << '\n';

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { { 0.0 }, { 5.0e-02 } };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );
    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime     = 0.10;
    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value() - tReferenceTime ) / tReferenceTime );
    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 1.334186887938226e+00;
    real tRelTempDifference    = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 2 ) - tReferenceTemperature ) / tReferenceTemperature );
    REQUIRE( tRelTempDifference < 1.0e-8 );

    // check x-displacement at node aNodeId in first time step
    real tReferenceXDisplacement = -8.310495618961509e-04;
    real tRelXDispDifference     = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 3, 2 ) - tReferenceXDisplacement ) / tReferenceXDisplacement );
    REQUIRE( tRelXDispDifference < 1.0e-6 );

    // check x-displacement at node aNodeId in first time step
    real tReferenceYDisplacement = 5.243129489358523e-03;
    real tRelYDispDifference     = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 4, 2 ) - tReferenceYDisplacement ) / tReferenceYDisplacement );
    REQUIRE( tRelYDispDifference < 1.0e-6 );

    // check volume-IQI at node aNodeId in first time step
    real tReferenceVolume     = 5.00e-05;
    real tRelVolumeDifference = std::abs( ( aExoIO.get_global_variable( 0, 0 ) - tReferenceVolume ) / tReferenceVolume );
    REQUIRE( tRelVolumeDifference < 1.0e-8 );
}

//---------------------------------------------------------------

extern "C" void
check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( "single_element.exo", 0, false, false );

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
        REQUIRE( tNumNodes == 4 );
        REQUIRE( tNumElems == 1 );
    }

    // check results
    uint tNodeId = 3;

    check_linear_results( tExoIO, tNodeId );
}

//---------------------------------------------------------------

extern "C" void
check_quadratic_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( "single_element.exo", 0, false, false );

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
        REQUIRE( tNumNodes == 9 );
        REQUIRE( tNumElems == 1 );
    }

    // check results
    uint tNodeId = 3;

    check_quadratic_results( tExoIO, tNodeId );
}

//---------------------------------------------------------------

TEST_CASE( "Thermo_Elastic_Element",
        "[moris],[example],[structure],[thermo_elastic],[Thermo_Elastic_Element],[single_element]" )
{
    // write out log file
    gLogger.initialize( "Log.log" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./single_element.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet == 0 );

    // set interpolation order
    gInterpolationOrder = 1;

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
