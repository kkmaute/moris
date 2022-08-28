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
void check_linear_results(moris::mtk::Exodus_IO_Helper & aExoIO,uint aNodeId)
{
    if (gPrintReferenceValues)
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point");

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision(15) << aExoIO.get_time_value() << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { {+7.600000000000002e-01},{+2.422727272727273e-01} };

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = aExoIO.get_nodal_coordinate( aNodeId );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);

    MORIS_LOG_INFO("Check nodal x-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate(0),tActualCoordinate(0),tRelDiffNorm*100.0);
    MORIS_LOG_INFO("Check nodal y-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate(1),tActualCoordinate(1),tRelDiffNorm*100.0);

    REQUIRE( tRelDiffNorm <  1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.000000000000000e+00;

    real tActualTime = aExoIO.get_time_value( );

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime) / tReferenceTime );

    MORIS_LOG_INFO("Check time:                 reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTime,tActualTime,tRelDiffNorm*100.0);

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 8.671607168734002e+04;

    real tActualTemperature = aExoIO.get_nodal_field_value( aNodeId, 2, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature ) / tReferenceTemperature );

    MORIS_LOG_INFO("Check nodal temperature:    reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTemperature,tActualTemperature,tRelTempDifference*100.0);

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelTempDifference < 1.0e-4);

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    real tReferenceIQI = 8.326070772043242e+04;

    real tActualIQI = aExoIO.get_global_variable(0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI ) / tReferenceIQI );

    MORIS_LOG_INFO("Check temperature IQI:      reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceIQI,tActualIQI,tRelIQIDifference*100.0);

    REQUIRE(  tRelIQIDifference < 1.0e-4);
}

//---------------------------------------------------------------

extern "C"
void check_quadratic_results(moris::mtk::Exodus_IO_Helper & aExoIO,uint aNodeId)
{
    if (gPrintReferenceValues)
    {
        // coordinates of reference point
        moris::print( aExoIO.get_nodal_coordinate( aNodeId ), "Coordinates of reference point");

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision(15) << aExoIO.get_time_value() << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { {+7.600000000000002e-01},{+2.422727272727273e-01} };

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = aExoIO.get_nodal_coordinate( aNodeId );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);

    MORIS_LOG_INFO("Check nodal x-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate(0),tActualCoordinate(0),tRelDiffNorm*100.0);
    MORIS_LOG_INFO("Check nodal y-coordinates:  reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceCoordinate(1),tActualCoordinate(1),tRelDiffNorm*100.0);

    REQUIRE( tRelDiffNorm <  1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.000000000000000e+00;

    real tActualTime = aExoIO.get_time_value( );

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime) / tReferenceTime );

    MORIS_LOG_INFO("Check time:                 reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTime,tActualTime,tRelDiffNorm*100.0);

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 8.673830875995508e+04;

    real tActualTemperature = aExoIO.get_nodal_field_value( aNodeId, 2, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature ) / tReferenceTemperature );

    MORIS_LOG_INFO("Check nodal temperature:    reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceTemperature,tActualTemperature,tRelTempDifference*100.0);

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelTempDifference < 1.0e-4);

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    real tReferenceIQI = 8.329556290215004e+04;

    real tActualIQI = aExoIO.get_global_variable(0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI ) / tReferenceIQI );

    MORIS_LOG_INFO("Check temperature IQI:      reference %12.5e, actual %12.5e, prec. error %12.5e.",
            tReferenceIQI,tActualIQI,tRelIQIDifference*100.0);

    REQUIRE(  tRelIQIDifference < 1.0e-4);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Linear - Serial.");
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static_Temp_Only.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, prec. error %12.5e.",2,tNumDims,std::abs((tNumDims-2)/2*100.0));
        MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, perc. error %12.5e.",6728,tNumNodes,std::abs((tNumNodes-5864)/5864*100.0));
        MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, perc. error %12.5e.",5754,tNumElems,std::abs((tNumElems-5754)/5754*100.0));

        REQUIRE( tNumDims  ==  2    );
        REQUIRE( tNumNodes ==  6728 );
        REQUIRE( tNumElems ==  5754 );
    }

    // check results
    uint tNodeId = 1798;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_parallel()
{
    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Linear - Parallel on %i processors.",par_size());
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static_Temp_Only.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, prec. error %12.5e.",2,tNumDims,std::abs((tNumDims-2)/2*100.0));
        MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, perc. error %12.5e.",4359,tNumNodes,std::abs((tNumNodes-3797)/3797*100.0));
        MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, perc. error %12.5e.",3714,tNumElems,std::abs((tNumElems-3714)/3714*100.0));

        REQUIRE( tNumDims  ==  2    );
        REQUIRE( tNumNodes ==  4359 );
        REQUIRE( tNumElems ==  3714 );
    }

    // check results
    uint tNodeId = 1117;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

extern "C"
void check_quadratic_results_serial()
{
    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Quadratic - Serial.");
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static_Temp_Only.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, prec. error %12.5e.",2,tNumDims,std::abs((tNumDims-2)/2*100.0));
        MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, perc. error %12.5e.",18711,tNumNodes,std::abs((tNumNodes-17847)/17847*100.0));
        MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, perc. error %12.5e.",6573,tNumElems,std::abs((tNumElems-6573)/6573*100.0));

        REQUIRE( tNumDims  ==  2     );
        REQUIRE( tNumNodes ==  18711 );
        REQUIRE( tNumElems ==  6573  );
    }

    // check results
    uint tNodeId = 8583;

    check_quadratic_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

extern "C"
void check_quadratic_results_parallel()
{
    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Checking Results - Quadratic - Parallel on %i processors.",par_size());
    MORIS_LOG_INFO("");

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static_Temp_Only.exo",0,false,false);

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    if (gPrintReferenceValues)
    {
        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
        std::cout << "Number of elements  : " << tNumElems << std::endl;
    }
    else
    {
        MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, prec. error %12.5e.",2,tNumDims,std::abs((tNumDims-2)/2*100.0));
        MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, perc. error %12.5e.",12065,tNumNodes,std::abs((tNumNodes-11503)/11503*100.0));
        MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, perc. error %12.5e.",4233,tNumElems,std::abs((tNumElems-4233)/4233*100.0));

        REQUIRE( tNumDims  ==  2     );
        REQUIRE( tNumNodes ==  12065 );
        REQUIRE( tNumElems ==  4233  );
    }

    // check results
    uint tNodeId = 5332;

    check_quadratic_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Channel_with_Four_Cylinders_Static_Linear",
        "[moris],[example],[thermal],[diffusion]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_with_Four_Cylinders_Static_Temp_Only.so";

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Executing Channel_with_Four_Cylinders_Static_Temp_Only: Interpolation order 1 - %i Processors.",par_size());
    MORIS_LOG_INFO("");

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            check_linear_results_serial();
            break;
        }
        case 2:
        {
            if (par_rank() == 1)
            {
                // set screen output processor
                gLogger.set_screen_output_rank(1);

                // perform check
                check_linear_results_parallel();

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

TEST_CASE("Channel_with_Four_Cylinders_Static_Quadratic",
        "[moris],[example],[thermal],[diffusion]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_with_Four_Cylinders_Static_Temp_Only.so";

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO("");
    MORIS_LOG_INFO("Executing Channel_with_Four_Cylinders_Static_Temp_Only: Interpolation order 2 - %i Processors.",par_size());
    MORIS_LOG_INFO("");

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    switch ( par_size() )
    {
        case 1:
        {
            check_quadratic_results_serial();
            break;
        }
        case 2:
        {
            if (par_rank() == 1)
            {
                // set screen output processor
                gLogger.set_screen_output_rank(1);

                // perform check
                check_quadratic_results_parallel();

                // set screen output processor
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

