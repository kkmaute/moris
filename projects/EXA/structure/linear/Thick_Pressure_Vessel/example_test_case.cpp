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
            "Pressure_Vessel_3D_Case_" + std::to_string( gCaseIndex ) + ".exo";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", gCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( tExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 11115, 8350, 11115, 3250, 3250 };

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
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 4, 0 ) << '\n';

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 5, 0 ) << '\n';

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 0, 0 ) << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 3, 3, 3, 3, 3 };
    Vector< uint > tReferenceNumNodes = { 17824, 6061, 17824, 5716, 5716 };
    Vector< uint > tReferenceNumElems = { 39168, 13034, 39168, 12288, 12288 };

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

    tReferenceCoordinate.push_back( { { +2.500000000000000e-01 }, { +2.500000000000000e-01 }, { +2.783105534096050e-01 } } );
    tReferenceCoordinate.push_back( { { +2.500000000000000e-01 }, { +2.500000000000000e-01 }, { +2.783105534096050e-01 } } );
    tReferenceCoordinate.push_back( { { +2.500000000000000e-01 }, { +2.500000000000000e-01 }, { +2.783105534096051e-01 } } );
    tReferenceCoordinate.push_back( { { +2.000000000000000e-01 }, { +2.000000000000000e-01 }, { +3.500000000000000e-01 } } );
    tReferenceCoordinate.push_back( { { +2.000000000000000e-01 }, { +2.000000000000000e-01 }, { +3.500000000000000e-01 } } );

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
    MORIS_LOG_INFO( "Check nodal z-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( gCaseIndex )( 2 ),
            tActualCoordinate( 2 ),
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelDiffNorm < 1.0e-5 );

    // check time value for time step index 0
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( gCaseIndex ) ) / tReferenceTime( gCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime( gCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceDisplacement;

    tReferenceDisplacement.push_back( { { -1.619383350258344e+00 }, { -1.619383350258375e+00 }, { -1.801581709307490e+00 } } );
    tReferenceDisplacement.push_back( { { -1.619383350258324e+00 }, { -1.619383350258371e+00 }, { -1.801581709307470e+00 } } );
    tReferenceDisplacement.push_back( { { -1.619383350258278e+00 }, { -1.619383350258378e+00 }, { -1.801581709307436e+00 } } );
    tReferenceDisplacement.push_back( { { -1.289499430891365e+00 }, { -1.289499430891467e+00 }, { -2.253455428719167e+00 } } );
    tReferenceDisplacement.push_back( { { -1.281913408757168e+00 }, { -1.281913408758300e+00 }, { -2.243085147963024e+00 } } );

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 4, 0 ) }
    };

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( gCaseIndex ) ) / norm( tReferenceDisplacement( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( gCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );

    // check temperature at node aNodeId in first time step (temperature is 6th nodal field, first time step has index 0)
    Vector< real > tReferenceTemperature;
    tReferenceTemperature.push_back( 2.000945133270288e+02 );
    tReferenceTemperature.push_back( 2.000945133270283e+02 );
    tReferenceTemperature.push_back( 2.000945133270288e+02 );
    tReferenceTemperature.push_back( 2.000221005039185e+02 );
    tReferenceTemperature.push_back( 2.003526463545086e+02 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 5, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature( gCaseIndex ) ) / tReferenceTemperature( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal temperature:    reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTemperature( gCaseIndex ),
            tActualTemperature,
            tRelTempDifference * 100.0 );

    REQUIRE( tRelTempDifference < 1.0e-5 );

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Vector< real > tReferenceIQI;
    tReferenceIQI.push_back( 4.911276929614884e-01 );
    tReferenceIQI.push_back( 4.911276929614933e-01 );
    tReferenceIQI.push_back( 4.911276929614822e-01 );
    tReferenceIQI.push_back( 4.911276929614835e-01 );
    tReferenceIQI.push_back( 4.916643196542578e-01 );

    real tActualIQI = tExoIO.get_global_variable( 0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI( gCaseIndex ) ) / tReferenceIQI( gCaseIndex ) );

    MORIS_LOG_INFO( "Check temperature IQI:      reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceIQI( gCaseIndex ),
            tActualIQI,
            tRelIQIDifference * 100.0 );

    REQUIRE( tRelIQIDifference < 1.0e-5 );
}

//---------------------------------------------------------------

TEST_CASE( "Pressure_Vessel_3D_Linear",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Pressure_Vessel_3D.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    // set case index
    gCaseIndex = par_size() == 1 ? 0 : 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Pressure_Vessel_3D: Interpolation order 1 - %i Processors.", par_size() );
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
            check_results();
            break;
        }
        // Test Case 1
        case 4:
        {
            if ( par_rank() == 1 )
            {
                // set screen output processor
                gLogger.set_screen_output_rank( 1 );

                // perform check
                check_results();

                // reset screen output processor
                gLogger.set_screen_output_rank( 0 );
            }
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Pressure_Vessel_3D_Immersed_Linear",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Pressure_Vessel_3D_Immersed.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    // set case index
    gCaseIndex = par_size() == 1 ? 2 : 3;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Pressure_Vessel_3D_Immersed: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 2
        case 1:
        {
            // perform check
            check_results();
            break;
        }
        // Test Case 3
        case 4:
        {
            if ( par_rank() == 1 )
            {
                // set screen output processor
                gLogger.set_screen_output_rank( 1 );

                // perform check
                check_results();

                // reset screen output processor
                gLogger.set_screen_output_rank( 0 );
            }
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Pressure_Vessel_3D_Immersed_Quadratic",
        "[moris],[example],[structure],[quadratic]" )
{
    // this test only runs in parallel; is skipped for serial runs
    if ( par_size() < 4 )
    {
        return;
    }

    // set case index
    gCaseIndex = 4;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Pressure_Vessel_3D_Immersed.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Pressure_Vessel_3D_Immersed: Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 4
        case 4:
        {
            if ( par_rank() == 1 )
            {
                // set screen output processor
                gLogger.set_screen_output_rank( 1 );

                // perform check
                check_results();

                // reset screen output processor
                gLogger.set_screen_output_rank( 0 );
            }
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}
