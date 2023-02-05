/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp.template
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

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        std::string aExoFileName,
        uint        aTestCaseIndex )
{
    MORIS_LOG_INFO( "" );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( "" );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, true, false );

    // define reference node IDs
    Cell< uint > tReferenceNodeId = { 729 };

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        // coordinates of reference point
        std::cout << "X coordinate: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) )( 0 ) << std::endl;
        std::cout << "Y coordinate: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) )( 1 ) << std::endl;

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << std::endl;

        // displacement values
        std::cout << "Displacement X-dir: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << std::endl;
        std::cout << "Displacement Y-dir: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 ) << std::endl;

        // stress at reference point
        std::cout << "Stress: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell< uint > tReferenceNumDims  = { 2 };
    Cell< uint > tReferenceNumNodes = { 7059 };
    Cell< uint > tReferenceNumElems = { 13228 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

    // define reference coordinates for node aNodeId
    Cell< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { 1.882561482489109e-02 }, { 2.055238932371140e-03 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( aTestCaseIndex ) ) / moris::norm( tReferenceCoordinate( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 0 ),
            tActualCoordinate( 0 ),
            tRelDiffNorm * 100.0 );
    MORIS_LOG_INFO( "Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 1 ),
            tActualCoordinate( 1 ),
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelDiffNorm < 1.0e-5 );

    // check time value for time step index 0
    Cell< real > tReferenceTime;
    tReferenceTime.push_back( 1 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check displacement at node aNodeId in first time step (displacements are 3rd and 4th nodal fields, first time step has index 0)
    Cell< Matrix< DDRMat > > tReferenceDisplacement;

    tReferenceDisplacement.push_back( { { -1.128321042350804e-04 }, { 9.853289535187406e-06 } } );

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) }
    };

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( aTestCaseIndex ) ) / norm( tReferenceDisplacement( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacement:   reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( aTestCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );

    // check temperature at node aNodeId in first time step (temperature is 7th nodal field, first time step has index 0)
    Cell< real > tReferenceTemperature;
    tReferenceTemperature.push_back( 1.038658447265625e+03 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature( aTestCaseIndex ) ) / tReferenceTemperature( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal temperature:    reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTemperature( aTestCaseIndex ),
            tActualTemperature,
            tRelTempDifference * 100.0 );

    REQUIRE( tRelTempDifference < 1.0e-5 );

    // check stresses at node aNodeId in first time step (stress is 5th & 6th nodal field)
    Cell< real > tReferenceStress;
    tReferenceStress.push_back( 2.397992760916115e+08 );

    real tActualStress = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 );

    real tRelStressDifference = std::abs( ( tActualStress - tReferenceStress( aTestCaseIndex ) ) / tReferenceStress( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal stress:         reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceStress( aTestCaseIndex ),
            tActualStress,
            tRelStressDifference * 100.0 );

    REQUIRE( tRelStressDifference < 1.0e-5 );
}

TEST_CASE( "AxisymmetricProblem",
        "[moris],[example],[structure],[thermo_elastic],[MOPAR]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./AxisymmetricProblem.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 0
        case 1:
        {
            // perform check
            check_results( "AxisymmetricProblem.exo", 0 );
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}
