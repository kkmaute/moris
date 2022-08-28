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
        std::cout << "X-Velocity at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) << std::endl;
        std::cout << "X-Velocity at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) << std::endl;
        std::cout << "Pressure at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) << std::endl;
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
                aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI - VX   value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(0, 0 ) << std::endl;
        std::cout << "IQI - VY   value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(1, 0 ) << std::endl;
        std::cout << "IQI - P    value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(2, 0 ) << std::endl;
        std::cout << "IQI - TEMP value: " << std::scientific << std::setprecision(15) << aExoIO.get_global_variable(3, 0 ) << std::endl;

        return;
    }

    // define reference coordinates for node aNodeId
    Matrix< DDRMat > tReferenceCoordinate = { {+7.600000000000002e-01},{+1.400000000000000e-01} };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);

    REQUIRE( tRelDiffNorm <  1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.000000000000000e+00;

    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value( ) - tReferenceTime) / tReferenceTime );

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceVelX = 0.5 * ( 1.593502502960416e+00 + 1.593502502960416e+00);
    real tReferenceVelY = 0.5 * (-3.831961099052318e-01 - 3.831961099052317e-01);
    real tReferencePres = 0.5 * ( 1.358249884899363e+00 + 1.358249884899363e+00);
    real tReferenceTemp = 0.5 * ( 1.371930816150207e+00 + 1.371930816150206e+00);

    real tRelDifference_VelX = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) - tReferenceVelX ) / tReferenceVelX );
    real tRelDifference_VelY = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 3, 0 ) - tReferenceVelY ) / tReferenceVelY );
    real tRelDifference_Pres = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 4, 0 ) - tReferencePres ) / tReferencePres );
    real tRelDifference_Temp = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 5, 0 ) - tReferenceTemp ) / tReferenceTemp );

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelDifference_VelX < 1.0e-5);
    REQUIRE(  tRelDifference_VelY < 1.0e-5);
    REQUIRE(  tRelDifference_Pres < 1.0e-5);
    REQUIRE(  tRelDifference_Temp < 1.0e-5);

    // check IQIs of first time step (only 1 IQI is defined, first time step has index 0)
    real tReferenceIQI_VelX = 0.5 * ( 9.015779490552553e-01 + 9.015779490552543e-01);
    real tReferenceIQI_VelY = 0.5 * (-1.581790313964251e-02 - 1.581790313964262e-02);
    real tReferenceIQI_Pres = 0.5 * ( 7.645059302522635e-01 + 7.645059302522605e-01);
    real tReferenceIQI_Temp = 0.5 * ( 1.372894021696889e+01 + 1.372894021696856e+01);

    real tRelIQIDifference_VelX = std::abs( ( aExoIO.get_global_variable(0, 0 ) - tReferenceIQI_VelX ) / tReferenceIQI_VelX );
    real tRelIQIDifference_VelY = std::abs( ( aExoIO.get_global_variable(1, 0 ) - tReferenceIQI_VelY ) / tReferenceIQI_VelY );
    real tRelIQIDifference_Pres = std::abs( ( aExoIO.get_global_variable(2, 0 ) - tReferenceIQI_Pres ) / tReferenceIQI_Pres );
    real tRelIQIDifference_Temp = std::abs( ( aExoIO.get_global_variable(3, 0 ) - tReferenceIQI_Temp ) / tReferenceIQI_Temp );

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelIQIDifference_VelX < 1.0e-5);
    REQUIRE(  tRelIQIDifference_VelY < 1.0e-5);
    REQUIRE(  tRelIQIDifference_Pres < 1.0e-5);
    REQUIRE(  tRelIQIDifference_Temp < 1.0e-5);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static.exo",0,false,false);

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

    REQUIRE( tNumDims  ==  2     );
    REQUIRE( tNumNodes ==  16887 );
    REQUIRE( tNumElems ==  15533 );

    // check results
    uint tNodeId = 3307;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_parallel()
{
    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO("Channel_with_Four_Cylinders_Static.exo",0,false,false);

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
        REQUIRE( tNumDims  ==  2    );
        REQUIRE( tNumNodes ==  3490 );
        REQUIRE( tNumElems ==  3244 );
    }

    // check results at reference node (watch: node Id depends on processor)
    uint tNodeId = 68;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

TEST_CASE("Channel_with_Four_Cylinders_Static_Linear",
        "[moris],[example],[thermal],[advection]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_with_Four_Cylinders_Static.so";

    char * argv[2] = {tString1,tString2};

    // set interpolation order
    gInterpolationOrder = 1;

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
        case 4:
        {
            if (par_rank() == 1)
            {
                check_linear_results_parallel();
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
        "[moris],[example],[thermal],[advection]")
{
    // quadratic case currently not working
    if (false)
    {
        // define command line call
        int argc = 2;

        char tString1[] = "";
        char tString2[] = "./Channel_with_Four_Cylinders_Static.so";

        char * argv[2] = {tString1,tString2};

        // set interpolation order
        gInterpolationOrder = 2;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // check results
        switch ( par_size() )
        {
            case 1:
            {
                // check_quadratic_results_serial();
                break;
            }
            case 4:
            {
                if (par_rank() == 1)
                {
                    // check_quadratic_results_parallel();
                }
                break;
            }
            default:
            {
                MORIS_ERROR(false,"Example problem not configured for %d processors.",par_size());
            }
        }
    }
}

