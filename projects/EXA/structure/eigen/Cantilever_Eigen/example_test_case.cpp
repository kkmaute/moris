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
uint gOrder;

// global variable for eigen-algorithm
std::string gPrecSolver;

// global test-case index
uint gTestCaseIndex;

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
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, true, false );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 14, 43 };

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
        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) ), "Coordinates of reference point" );

        // Displacement in X-direction at reference point
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << std::endl;

        // Displacement in Y-direction at reference point
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2 };
    Vector< uint > tReferenceNumNodes = { 49241, 49241 };
    Vector< uint > tReferenceNumElems = { 48000, 48000 };

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
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { 0.00175 }, { 0.00000 } } );
    tReferenceCoordinate.push_back( { { 0.00525 }, { 0.00025 } } );

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

    // check temperature at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceDisplacement;

    if ( aTestCaseIndex == 0 )
    {
        tReferenceDisplacement.push_back( { { -1.42942084733364e-06 }, { 6.46851986039073e-08 } } );
        tReferenceDisplacement.push_back( { { 0.000993543881481338 }, { -0.000532691066062199 } } );
    }
    else if ( aTestCaseIndex == 1 )
    {
        tReferenceDisplacement.push_back( { { -1.42942084733364e-06 }, { 6.46851986039073e-08 } } );
        tReferenceDisplacement.push_back( { { 0.000993543881481338 }, { -0.000532691066062199 } } );
    }

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) }
    };

    // adjust sign as computed eigen vector might be pointing in opposite direction than reference vector
    if ( tActualDisplacement( 0 ) * tReferenceDisplacement( aTestCaseIndex )( 0 ) < 0 )
    {
        tActualDisplacement = -1.0 * tActualDisplacement;
    }

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( aTestCaseIndex ) ) / norm( tReferenceDisplacement( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacement:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( aTestCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Cantilever_Eigen_Pardiso",
        "[moris],[example],[structure],[Pardiso]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Cantilever_Eigen.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing EigenProblem: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set eigen algorithm
    gPrecSolver = "Amesos_Pardiso";

    // set test-case index
    gTestCaseIndex = 0;

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    check_results( "Cantilever_Eigen_Pardiso.exo", gTestCaseIndex );
}

//------------------------------------------------------------------------------------------------------------------------------------------------

TEST_CASE( "Cantilever_Eigen_Umfpack",
        "[moris],[example],[structure],[Umfpack]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Cantilever_Eigen.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing EigenProblem: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set eigen algorithm
    gPrecSolver = "Amesos_Umfpack";

    // set test-case index
    gTestCaseIndex = 1;

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check test statement should follow
    REQUIRE( tRet == 0 );

    // Perform check results for test-case 0
    check_results( "Cantilever_Eigen_Umfpack.exo", gTestCaseIndex );
}
