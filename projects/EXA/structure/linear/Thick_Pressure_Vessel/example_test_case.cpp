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

#include "cl_Logger.hpp"                // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------

// global variable for interpolation order
uint gInterpolationOrder;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

extern "C"
void check_results(
        std::string aExoFileName,
        uint        aTestCaseIndex)
{
    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Test Case %d on %i processor.",aTestCaseIndex,par_size());
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO(aExoFileName.c_str(),0,false,false);

    // define reference node IDs
    Cell<uint> tReferenceNodeId  = {11319,4070,1658,584,584};

    if (gPrintReferenceValues)
    {
        std::cout << "Test case index: " << aTestCaseIndex << std::endl;

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;

        // coordinates of reference point
        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId(aTestCaseIndex) ), "Coordinates of reference point");

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision(15) << tExoIO.get_time_value() << std::endl;

        // solution of reference point at reference time step
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision(15) <<
                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 2, 0 ) << "," <<
                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 3, 0 ) << "," <<
                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) <<
                std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 5, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision(15) << tExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell<uint> tReferenceNumDims  = {3,     3,     3,     3,     3    };
    Cell<uint> tReferenceNumNodes = {18292, 6217,  17824, 5716,  5716 };
    Cell<uint> tReferenceNumElems = {40254, 13396, 39168, 12288, 12288};

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims(aTestCaseIndex),tNumDims,std::abs((tNumDims-tReferenceNumDims(aTestCaseIndex))/tReferenceNumDims(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes(aTestCaseIndex),tNumNodes,std::abs((tNumNodes-tReferenceNumNodes(aTestCaseIndex))/tReferenceNumNodes(aTestCaseIndex)*100.0));
    MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems(aTestCaseIndex),tNumElems,std::abs((tNumElems-tReferenceNumElems(aTestCaseIndex))/tReferenceNumElems(aTestCaseIndex)*100.0));

    REQUIRE( tNumDims  ==  tReferenceNumDims(aTestCaseIndex)  );
    REQUIRE( tNumNodes ==  tReferenceNumNodes(aTestCaseIndex) );
    REQUIRE( tNumElems ==  tReferenceNumElems(aTestCaseIndex) );

    // define reference coordinates for node aNodeId
    Cell<Matrix< DDRMat >> tReferenceCoordinate;

    tReferenceCoordinate.push_back( { {+2.500000000000000e-01},{ +2.500000000000000e-01},{+2.783105534096050e-01} } );
    tReferenceCoordinate.push_back( { {+2.500000000000000e-01},{ +2.500000000000000e-01},{+2.783105534096050e-01} } );
    tReferenceCoordinate.push_back( { {+2.000000000000000e-01},{ +2.000000000000000e-01},{+3.500000000000000e-01} } );
    tReferenceCoordinate.push_back( { {+2.000000000000000e-01},{ +2.000000000000000e-01},{+3.500000000000000e-01} } );
    tReferenceCoordinate.push_back( { {+2.000000000000000e-01},{ +2.000000000000000e-01},{+3.500000000000000e-01} } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId(aTestCaseIndex) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate(aTestCaseIndex) )/ moris::norm(tReferenceCoordinate(aTestCaseIndex));

    MORIS_LOG_INFO("Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate(aTestCaseIndex)(0),tActualCoordinate(0),tRelDiffNorm*100.0);
    MORIS_LOG_INFO("Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate(aTestCaseIndex)(1),tActualCoordinate(1),tRelDiffNorm*100.0);
    MORIS_LOG_INFO("Check nodal z-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate(aTestCaseIndex)(2),tActualCoordinate(2),tRelDiffNorm*100.0);

    REQUIRE( tRelDiffNorm <  1.0e-5 );

    // check time value for time step index 0
    Cell<real> tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );
    tReferenceTime.push_back( 1.000000000000000e+00 );

    real tActualTime = tExoIO.get_time_value( );

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime(aTestCaseIndex)) / tReferenceTime(aTestCaseIndex) );

    MORIS_LOG_INFO("Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime(aTestCaseIndex),tActualTime,tRelDiffNorm*100.0);

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
     Cell<Matrix< DDRMat >> tReferenceDisplacement;

     tReferenceDisplacement.push_back( { {-1.619414637604386e+00},{-1.619414604667233e+00},{-1.801613073132456e+00} } );
     tReferenceDisplacement.push_back( { {-1.619414637604319e+00},{-1.619414604667178e+00},{-1.801613073132412e+00} } );
     tReferenceDisplacement.push_back( { {-1.289499430891421e+00},{-1.289499430891530e+00},{-2.253455428719179e+00} } );
     tReferenceDisplacement.push_back( { {-1.289499430891418e+00},{-1.289499430891522e+00},{-2.253455428719149e+00} } );
     tReferenceDisplacement.push_back( { {-1.281913408756792e+00},{-1.281913408758628e+00},{-2.243085147962265e+00} } );

     Matrix< DDRMat > tActualDisplacement = {
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 2, 0 ) },
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 3, 0 ) },
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) } };

     real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement(aTestCaseIndex) ) / norm( tReferenceDisplacement(aTestCaseIndex) );

     MORIS_LOG_INFO("Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
             norm(tReferenceDisplacement(aTestCaseIndex)) ,norm(tActualDisplacement),tRelDispDifference*100.0);

    REQUIRE(  tRelDispDifference < 1.0e-5);

    // check temperature at node aNodeId in first time step (temperature is 6th nodal field, first time step has index 0)
    Cell<real> tReferenceTemperature;
    tReferenceTemperature.push_back( 2.000944820689973e+02 );
    tReferenceTemperature.push_back( 2.000944820689973e+02 );
    tReferenceTemperature.push_back( 2.000221005039182e+02 );
    tReferenceTemperature.push_back( 2.000221005039181e+02 );
    tReferenceTemperature.push_back( 2.003526463545087e+02 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 5, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature(aTestCaseIndex) ) / tReferenceTemperature(aTestCaseIndex) );

    MORIS_LOG_INFO("Check nodal temperature:    reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTemperature(aTestCaseIndex),tActualTemperature,tRelTempDifference*100.0);

    REQUIRE(  tRelTempDifference < 1.0e-5);

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Cell<real> tReferenceIQI;
    tReferenceIQI.push_back( 4.911266018484971e-01 );
    tReferenceIQI.push_back( 4.911266018484888e-01 );
    tReferenceIQI.push_back( 4.911276929614986e-01 );
    tReferenceIQI.push_back( 4.911276929614983e-01 );
    tReferenceIQI.push_back( 4.916643196541598e-01 );

    real tActualIQI = tExoIO.get_global_variable(0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI(aTestCaseIndex) ) / tReferenceIQI(aTestCaseIndex) );

    MORIS_LOG_INFO("Check temperature IQI:      reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceIQI(aTestCaseIndex),tActualIQI,tRelIQIDifference*100.0);

    REQUIRE(  tRelIQIDifference < 1.0e-5);
}

//---------------------------------------------------------------

TEST_CASE("Pressure_Vessel_3D_Linear",
        "[moris],[example],[structure],[linear]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Pressure_Vessel_3D.so";

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Executing Pressure_Vessel_3D: Interpolation order 1 - %i Processors.",par_size());
    MORIS_LOG_INFO("");

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 0
        case 1:
        {
            // perform check
            check_results("Pressure_Vessel_3D.exo",0);
            break;
        }
        // Test Case 1
        case 4:
        {
            if (par_rank() == 1)
            {
                // set screen output processor
                gLogger.set_screen_output_rank(1);

                // perform check
                check_results("Pressure_Vessel_3D.exo",1);

                // reset screen output processor
                gLogger.set_screen_output_rank(0);
            }
            break;
        }
        default:
        {
            MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
        }
    }
}

//---------------------------------------------------------------

TEST_CASE("Pressure_Vessel_3D_Immersed_Linear",
        "[moris],[example],[structure],[linear]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Pressure_Vessel_3D_Immersed.so";

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Executing Pressure_Vessel_3D_Immersed: Interpolation order 1 - %i Processors.",par_size());
    MORIS_LOG_INFO("");

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 2
        case 1:
        {
            // perform check
            check_results("Pressure_Vessel_3D_Immersed.exo",2);
            break;
        }
        // Test Case 3
        case 4:
        {
            if (par_rank() == 1)
            {
                // set screen output processor
                gLogger.set_screen_output_rank(1);

                // perform check
                check_results("Pressure_Vessel_3D_Immersed.exo",3);

                // reset screen output processor
                gLogger.set_screen_output_rank(0);
            }
            break;
        }
        default:
        {
            MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
        }
    }
}

//---------------------------------------------------------------

TEST_CASE("Pressure_Vessel_3D_Immersed_Quadratic",
        "[moris],[example],[structure],[quadratic]")
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

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Executing Pressure_Vessel_3D_Immersed: Interpolation order 2 - %i Processors.",par_size());
    MORIS_LOG_INFO("");

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        // Test Case 4
        case 4:
        {
            if (par_rank() == 1)
            {
                // set screen output processor
                gLogger.set_screen_output_rank(1);

                // perform check
                check_results("Pressure_Vessel_3D_Immersed.exo",4);

                // reset screen output processor
                gLogger.set_screen_output_rank(0);
            }
            break;
        }
        default:
        {
            MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
        }
    }
}

