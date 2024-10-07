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
check_results( moris::mtk::Exodus_IO_Helper &aExoIO, uint aNodeId )
{
    if ( gPrintReferenceValues )
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "X-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << '\n';
        std::cout << "Y-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << '\n';
        std::cout << "Z-Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << '\n';
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) << '\n';

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { { +2.000000000000000e00 }, { +1.250000000000000e-01 }, { +1.250000000000000e-01 } };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate ) / moris::norm( tReferenceCoordinate );

    REQUIRE( tRelDiffNorm < 1.0e-12 );

    // check time value for time step index 0
    real tReferenceTime = 1.000000000000000e+00;

    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value() - tReferenceTime ) / tReferenceTime );

    REQUIRE( tRelTimeDifference < 1.0e-12 );

    // check nodal val
    real tReferenceUX   = 2.003084205519877e+00;
    real tReferenceUY   = -4.907031775787233e-13;
    real tReferenceUZ   = 1.806606057303358e-12;
    real tReferenceTemp = 3.000000000000023e+00;

    real tRelDifference_UX   = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) - tReferenceUX ) / tReferenceUX );
    real tRelDifference_UY   = std::abs( aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) - tReferenceUY );
    real tRelDifference_UZ   = std::abs( aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) - tReferenceUZ );
    real tRelDifference_Temp = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) - tReferenceTemp ) / tReferenceTemp );

    // FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE( tRelDifference_UX < 1.0e-6 );
    REQUIRE( tRelDifference_UY < 1.0e-8 );
    REQUIRE( tRelDifference_UZ < 1.0e-8 );
    REQUIRE( tRelDifference_Temp < 1.0e-6 );
}

//---------------------------------------------------------------

extern "C" void
check_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( "Bar_Lin_Disp_Quad_Diff.exo", 0, false, false );

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

    REQUIRE( tNumDims == 3 );
    REQUIRE( tNumNodes == 390 );
    REQUIRE( tNumElems == 332 );

    // check results
    uint tNodeId = 232;

    check_results( tExoIO, tNodeId );
}

//---------------------------------------------------------------

TEST_CASE( "Thermoelastic_Bar_Linear_Disp_Quadratic_Diff",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Bar_Lin_Disp_Quad_Diff.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Bar2D: Displacement Interpolation order 1 - Diffusion Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 0
        case 1:
        {
            // perform check
            check_results_serial();
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

//---------------------------------------------------------------
