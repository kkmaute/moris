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
check_results(
        std::string aExoFileName,
        uint        aTestCaseIndex )
{

    MORIS_LOG_INFO( "" );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( "" );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        //        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        //        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        //        std::cout << "Number of elements  : " << tNumElems << std::endl;
        //
        //        // coordinates of reference point
        //        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point");
        //
        //        // time value for reference time step
        //        std::cout << "Time value: " << std::scientific << std::setprecision(15) << aExoIO.get_time_value() << std::endl;
        //
        //        // solution of reference point at reference time step
        //        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
        //                aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        //
        //        // value of IQI at reference time step
        //        std::cout << "IQI value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell< uint > tReferenceNumDims  = { 3, 3, 3, 3 };
    Cell< uint > tReferenceNumNodes = { 0, 0, 0, 0 };
    Cell< uint > tReferenceNumElems = { 0, 0, 0, 0 };

    Cell< uint > tReferenceNodeId = { 0, 0, 0, 0 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, prec. error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, perc. error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, perc. error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - ( aTestCaseIndex ) ) / (aTestCaseIndex)*100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumNodes( aTestCaseIndex ) );

    // define reference coordinates for node aNodeId
    Cell< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { +7.600000000000002e-01 }, { +2.422727272727273e-01 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( aTestCaseIndex ) ) / moris::norm( tReferenceCoordinate( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal x-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 0 ),
            tActualCoordinate( 0 ),
            tRelDiffNorm * 100.0 );
    MORIS_LOG_INFO( "Check nodal y-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 1 ),
            tActualCoordinate( 1 ),
            tRelDiffNorm * 100.0 );
    MORIS_LOG_INFO( "Check nodal z-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate( aTestCaseIndex )( 2 ),
            tActualCoordinate( 2 ),
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    Cell< real > tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+00 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    Cell< real > tReferenceTemperature;
    tReferenceTemperature.push_back( 8.670906164633136e+04 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature( aTestCaseIndex ) ) / tReferenceTemperature( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal temperature:    reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTemperature( aTestCaseIndex ),
            tActualTemperature,
            tRelTempDifference * 100.0 );

    // FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE( tRelTempDifference < 1.0e-4 );

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Cell< real > tReferenceIQI;
    tReferenceIQI.push_back( 8.324886510380027e+04 );

    real tActualIQI = tExoIO.get_global_variable( 0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI( aTestCaseIndex ) ) / tReferenceIQI( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check temperature IQI:      reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceIQI( aTestCaseIndex ),
            tActualIQI,
            tRelIQIDifference * 100.0 );

    REQUIRE( tRelIQIDifference < 1.0e-4 );
}

//---------------------------------------------------------------

TEST_CASE( "Heated_Sphere_Linear",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Heated_Sphere.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( "" );
    MORIS_LOG_INFO( "Executing Heated_Sphere: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( "" );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 0
        case 1:
        {
            // perform check
            // check_results("Heated_Sphere.exo",0);
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
                // check_results("Heated_Sphere.exo.1",1);

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

// TEST_CASE("Heated_Sphere_Quadratic",
//         "[moris],[example],[structure],[quadratic]")
//{
//     // define command line call
//     int argc = 2;
//
//     char tString1[] = "";
//     char tString2[] = "./Heated_Sphere.so";
//
//     char * argv[2] = {tString1,tString2};
//
//     // set interpolation order
//     gInterpolationOrder = 2;
//
//     MORIS_LOG_INFO("");
//     MORIS_LOG_INFO("Executing Heated_Sphere: Interpolation order 2 - %i Processors.",par_size());
//     MORIS_LOG_INFO("");
//
//     // call to performance manager main interface
//     fn_WRK_Workflow_Main_Interface( argc, argv );
//
//     // check results
//     switch ( par_size() )
//     {
//         // Test Case 2
//         case 1:
//         {
//             // perform check
//             check_results("Heated_Sphere.exo",2);
//             break;
//         }
//         // Test Case 3
//         case 4:
//         {
//             if (par_rank() == 1)
//             {
//                 // set screen output processor
//                 gLogger.set_screen_output_rank(1);
//
//                 // perform check
//                 //check_results("Heated_Sphere.exo.1",3);
//
//                 // reset screen output processor
//                 gLogger.set_screen_output_rank(0);
//             }
//             break;
//         }
//         default:
//         {
//             MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
//         }
//     }
// }
