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
#include "fn_stringify_matrix.hpp"

using namespace moris;

//---------------------------------------------------------------

// global variable for interpolation order
uint gInterpolationOrder;

// problem dimension: 2D or 3D
uint gDim;

// test case index
uint gTestCaseIndex;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

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
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 41, 72, 66, 129 };

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

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI 0 value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 0, 0 ) << std::endl;
        std::cout << "IQI 1 value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 1, 0 ) << std::endl;
        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 3, 2, 3 };
    Vector< uint > tReferenceNumNodes = { 45, 235, 70, 390 };
    Vector< uint > tReferenceNumElems = { 26, 332, 26, 332 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

    // define reference coordinates for node aNodeId
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { +4.1 }, { +0.5 } } );
    tReferenceCoordinate.push_back( { { +4.1 }, { +0.5 }, { +0.5 } } );
    tReferenceCoordinate.push_back( { { +4.1 }, { +0.5 } } );
    tReferenceCoordinate.push_back( { { +4.1 }, { +0.5 }, { +0.5 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( aTestCaseIndex ) ) / moris::norm( tReferenceCoordinate( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 0 ),
            tActualCoordinate( 0 ),
            tRelDiffNorm * 100.0 );
    MORIS_LOG_INFO( "Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 1 ),
            tActualCoordinate( 1 ),
            tRelDiffNorm * 100.0 );

    if ( tNumDims == 3 )
    {
        MORIS_LOG_INFO( "Check nodal z-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
                tReferenceCoordinate( aTestCaseIndex )( 2 ),
                tActualCoordinate( 2 ),
                tRelDiffNorm * 100.0 );
    }

    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    Vector< real > tReferenceTemperature;
    tReferenceTemperature.push_back( 7.731885467530020e+01 );
    tReferenceTemperature.push_back( 7.729432195240925e+01 );
    tReferenceTemperature.push_back( 7.728500000000312e+01 );
    tReferenceTemperature.push_back( 7.728500000001044e+01 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature( aTestCaseIndex ) ) / tReferenceTemperature( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal temperature:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTemperature( aTestCaseIndex ),
            tActualTemperature,
            tRelTempDifference * 100.0 );

    // FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE( tRelTempDifference < 1.0e-4 );

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Vector< real > tReferenceIQI_0;
    tReferenceIQI_0.push_back( 3.481332925575233e+02 );
    tReferenceIQI_0.push_back( 3.479932121344261e+02 );
    tReferenceIQI_0.push_back( 8.922568376614086e-22 );
    tReferenceIQI_0.push_back( 5.762418056802439e-20 );

    Vector< real > tReferenceIQI_1;
    tReferenceIQI_1.push_back( 6.838626208796658e+02 );
    tReferenceIQI_1.push_back( 6.838626251218382e+02 );
    tReferenceIQI_1.push_back( 2.911872014125403e-23 );
    tReferenceIQI_1.push_back( 2.805634603810048e-21 );

    real tActualIQI_0 = tExoIO.get_global_variable( 0, 0 );
    real tActualIQI_1 = tExoIO.get_global_variable( 1, 0 );

    real tAbsIQIDifference_0 = std::abs( ( tActualIQI_0 - tReferenceIQI_0( aTestCaseIndex ) ) );
    real tAbsIQIDifference_1 = std::abs( ( tActualIQI_1 - tReferenceIQI_1( aTestCaseIndex ) ) );

    MORIS_LOG_INFO( "Check temperature IQI 0:    reference %12.5e, actual %12.5e, absolute error %12.5e.",
            tReferenceIQI_0( aTestCaseIndex ),
            tActualIQI_0,
            tAbsIQIDifference_0 );
    MORIS_LOG_INFO( "Check temperature IQI 1:    reference %12.5e, actual %12.5e, absolute error %12.5e.",
            tReferenceIQI_1( aTestCaseIndex ),
            tActualIQI_1,
            tAbsIQIDifference_1 );

    REQUIRE( tAbsIQIDifference_0 < 1.0e-6 );
    REQUIRE( tAbsIQIDifference_1 < 1.0e-6 );
}

//---------------------------------------------------------------

extern "C" void
check_results_xtk_mesh(
        std::string aExoFileNameIG,
        std::string aExoFileNameIP )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case on %i processor.", par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIOIGMesh( aExoFileNameIG.c_str(), 0, false, false );
    moris::mtk::Exodus_IO_Helper tExoIOIPMesh( aExoFileNameIP.c_str(), 0, false, false );

    // define reference node IDs ( first is the IG mesh, second is the IP mesh )
    Vector< uint > tReferenceNodeId = { 8, 5 };

    if ( gPrintReferenceValues )
    {
        uint tNumDimsIG  = tExoIOIGMesh.get_number_of_dimensions();
        uint tNumNodesIG = tExoIOIGMesh.get_number_of_nodes();
        uint tNumElemsIG = tExoIOIGMesh.get_number_of_elements();

        uint tNumDimsIP  = tExoIOIPMesh.get_number_of_dimensions();
        uint tNumNodesIP = tExoIOIPMesh.get_number_of_nodes();
        uint tNumElemsIP = tExoIOIPMesh.get_number_of_elements();

        std::cout << "Number of dimensions IG mesh : " << tNumDimsIG << std::endl;
        std::cout << "Number of nodes IG mesh      : " << tNumNodesIG << std::endl;
        std::cout << "Number of elements IG mesh   : " << tNumElemsIG << std::endl;

        std::cout << "Number of dimensions IP mesh : " << tNumDimsIP << std::endl;
        std::cout << "Number of nodes IP mesh      : " << tNumNodesIP << std::endl;
        std::cout << "Number of elements IP mesh   : " << tNumElemsIP << std::endl;

        // coordinates of reference point
        moris::print( tExoIOIGMesh.get_nodal_coordinate( tReferenceNodeId( 0 ) ), "Coordinates of reference point IG mesh" );
        moris::print( tExoIOIPMesh.get_nodal_coordinate( tReferenceNodeId( 1 ) ), "Coordinates of reference point IP mesh" );

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2 };
    Vector< uint > tReferenceNumNodes = { 54, 60 };
    Vector< uint > tReferenceNumElems = { 26, 10 };

    uint tNumDimsIG  = tExoIOIGMesh.get_number_of_dimensions();
    uint tNumNodesIG = tExoIOIGMesh.get_number_of_nodes();
    uint tNumElemsIG = tExoIOIGMesh.get_number_of_elements();

    uint tNumDimsIP  = tExoIOIPMesh.get_number_of_dimensions();
    uint tNumNodesIP = tExoIOIPMesh.get_number_of_nodes();
    uint tNumElemsIP = tExoIOIPMesh.get_number_of_elements();

    REQUIRE( tNumDimsIG == tReferenceNumDims( 0 ) );
    REQUIRE( tNumNodesIG == tReferenceNumNodes( 0 ) );
    REQUIRE( tNumElemsIG == tReferenceNumElems( 0 ) );

    REQUIRE( tNumDimsIP == tReferenceNumDims( 1 ) );
    REQUIRE( tNumNodesIP == tReferenceNumNodes( 1 ) );
    REQUIRE( tNumElemsIP == tReferenceNumElems( 1 ) );

    // define reference coordinates for node aNodeId
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { +1.25 }, { +0.25 } } );
    tReferenceCoordinate.push_back( { { +1.25 }, { +0.0 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinateIG = tExoIOIGMesh.get_nodal_coordinate( tReferenceNodeId( 0 ) );
    Matrix< DDRMat > tActualCoordinateIP = tExoIOIPMesh.get_nodal_coordinate( tReferenceNodeId( 1 ) );

    real tRelDiffNormIG = moris::norm( tActualCoordinateIG - tReferenceCoordinate( 0 ) ) / moris::norm( tReferenceCoordinate( 0 ) );
    real tRelDiffNormIP = moris::norm( tActualCoordinateIP - tReferenceCoordinate( 1 ) ) / moris::norm( tReferenceCoordinate( 1 ) );

    MORIS_LOG_INFO( "Check nodal IG-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
            norm( tReferenceCoordinate( 0 ) ),
            norm( tActualCoordinateIG ),
            tRelDiffNormIG * 100.0 );
    MORIS_LOG_INFO( "Check nodal IP-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
            norm( tReferenceCoordinate( 1 ) ),
            norm( tActualCoordinateIP ),
            tRelDiffNormIP * 100.0 );

    // Define reference values for nodal fields
    Matrix< DDRMat > tReferenceValuesIGMesh = { { 1.0, -10.0 } };
    Matrix< DDRMat > tReferenceValuesIPMesh = { { 0.0625, 0.375, 0.0625 } };

    // Define actual values for nodal fields
    Matrix< DDRMat > tActualValuesIGMesh( 1, 2, 0.0 );
    tActualValuesIGMesh( 0 ) = tExoIOIGMesh.get_nodal_field_value( tReferenceNodeId( 0 ), "weights_mi_0_ind_4", 0.0 );
    tActualValuesIGMesh( 1 ) = tExoIOIGMesh.get_nodal_field_value( tReferenceNodeId( 0 ), "weights_mi_0_ind_22", 0.0 );

    Matrix< DDRMat > tActualValuesIPMesh( 1, 3, 0.0 );
    tActualValuesIPMesh( 0 ) = tExoIOIPMesh.get_nodal_field_value( tReferenceNodeId( 1 ), "Basis_mi_0_ind_0", 0.0 );
    tActualValuesIPMesh( 1 ) = tExoIOIPMesh.get_nodal_field_value( tReferenceNodeId( 1 ), "Basis_mi_0_ind_4", 0.0 );
    tActualValuesIPMesh( 2 ) = tExoIOIPMesh.get_nodal_field_value( tReferenceNodeId( 1 ), "weights_mi_0_ind_10", 0.0 );

    // define lambdas to check if values are equal
    bool tIGMeshValuesAreEqual = std::equal( tActualValuesIGMesh.begin(), tActualValuesIGMesh.end(), tReferenceValuesIGMesh.begin() );
    bool tIPMeshValuesAreEqual = std::equal( tActualValuesIPMesh.begin(), tActualValuesIPMesh.end(), tReferenceValuesIPMesh.begin() );

    // print out IP and IG mesh values
    MORIS_LOG_INFO( "Check IG mesh values: reference %s , actual %s.", ios::stringify_log( tReferenceValuesIGMesh ).c_str(), ios::stringify_log( tActualValuesIGMesh ).c_str() );
    MORIS_LOG_INFO( "Check IP mesh values: reference %s , actual %s.", ios::stringify_log( tReferenceValuesIPMesh ).c_str(), ios::stringify_log( tActualValuesIPMesh ).c_str() );

    // check the values are matching
    REQUIRE( tIGMeshValuesAreEqual );
    REQUIRE( tIPMeshValuesAreEqual );
}

//---------------------------------------------------------------

TEST_CASE( "Two_Material_Bar_Linear",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Two_Material_Bar.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Two_Material_Bar - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 0;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 0
    check_results( "Two_Material_Bar_0.exo", gTestCaseIndex );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Two_Material_Bar - 3D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 3D
    gDim = 3;

    // set test case index
    gTestCaseIndex = 1;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 1
    check_results( "Two_Material_Bar_1.exo", gTestCaseIndex );
}

//---------------------------------------------------------------

TEST_CASE( "Two_Material_Bar_Quadratic",
        "[moris],[example],[structure],[quadratic]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Two_Material_Bar.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Two_Material_Bar - 2D: Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 2;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 2
    check_results( "Two_Material_Bar_2.exo", gTestCaseIndex );

#ifdef MORIS_HAVE_DEBUG
    // check xtk mesh results
    check_results_xtk_mesh( "xtk_temp.exo", "xtk_temp_ip.exo" );
#endif

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Two_Material_Bar - 3D: Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 3D
    gDim = 3;

    // set test case index
    gTestCaseIndex = 3;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 3
    check_results( "Two_Material_Bar_3.exo", gTestCaseIndex );
}
