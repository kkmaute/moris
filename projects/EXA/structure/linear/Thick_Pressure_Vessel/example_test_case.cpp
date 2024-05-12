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
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 11319, 8550, 1658, 3250, 3250 };

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << std::endl;

        // solution of reference point at reference time step
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 5, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 0, 0 ) << std::endl;

        return;
    }

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Vector< real > tReferenceIQI;
    tReferenceIQI.push_back( 4.911266018484971e-01 );
    tReferenceIQI.push_back( 4.911266018484888e-01 );
    tReferenceIQI.push_back( 4.911276929614986e-01 );
    tReferenceIQI.push_back( 4.911276929614983e-01 );
    tReferenceIQI.push_back( 4.916643196541598e-01 );

    real tActualIQI = tExoIO.get_global_variable( 0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI( aTestCaseIndex ) ) / tReferenceIQI( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check temperature IQI:      reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceIQI( aTestCaseIndex ),
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
            check_results( "Pressure_Vessel_3D.exo", 0 );
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
                check_results( "Pressure_Vessel_3D.exo", 1 );

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
            check_results( "Pressure_Vessel_3D_Immersed.exo", 2 );
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
                check_results( "Pressure_Vessel_3D_Immersed.exo", 3 );

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
                check_results( "Pressure_Vessel_3D_Immersed.exo", 4 );

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
