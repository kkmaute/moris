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

// text case index
uint gCaseIndex;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------

extern "C" void check_results()
{
    std::string tExoFileName =
            "Parabolic_Indenter_Linear_Case_" + std::to_string( gCaseIndex ) + ".e-s.0000";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", gCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( tExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 294 };

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << gCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        // coordinates of reference point
        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId( gCaseIndex ) ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2 };
    Vector< uint > tReferenceNumNodes = { 313 };
    Vector< uint > tReferenceNumElems = { 192 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims( gCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( gCaseIndex ) ) / tReferenceNumDims( gCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes( gCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( gCaseIndex ) ) / tReferenceNumNodes( gCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems( gCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( gCaseIndex ) ) / tReferenceNumElems( gCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( gCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( gCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( gCaseIndex ) );

    // define reference coordinates for node aNodeId
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { +5.000000000000001e-01 }, { +5.025000000000001e-01 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( gCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( gCaseIndex ) ) / moris::norm( tReferenceCoordinate( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( gCaseIndex )( 0 ),
            tActualCoordinate( 0 ),
            tRelDiffNorm * 100.0 );
    MORIS_LOG_INFO( "Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( gCaseIndex )( 1 ),
            tActualCoordinate( 1 ),
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelDiffNorm < 1.0e-5 );

    // check time value for time step index 0
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+01 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( gCaseIndex ) ) / tReferenceTime( gCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime( gCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceDisplacement;

    tReferenceDisplacement.push_back( { { -3.534434707883998e-06 }, { -2.785067809916309e-02 } } );

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) },
    };

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( gCaseIndex ) ) / norm( tReferenceDisplacement( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( gCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );
}

//---------------------------------------------------------------

TEST_CASE( "Parabolic_Indenter_Linear",
        "[moris],[example],[structure],[linear]" )
{
#ifdef MORIS_HAVE_ARBORX
    // check that run is serial; parallel not implemented yet
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Parabolic_Indenter_Linear.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    // set case index
    gCaseIndex = 0;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Parabolic_Indenter_Linear: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    check_results();

#else
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO("Parabolic_Indenter_Linear: Example skipped as Arborx not installed");
    MORIS_LOG_INFO( " " );
#endif
}
