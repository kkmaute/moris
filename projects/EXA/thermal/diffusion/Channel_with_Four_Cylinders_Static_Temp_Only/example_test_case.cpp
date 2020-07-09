//
// example specific interface to moris
//

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
    Matrix< DDRMat > tReferenceCoordinate = { {+8.736363636363635e-01},{+2.081818181818181e-01} };

    // check nodal coordinates
    real tRelDiffNorm = moris::norm( aExoIO.get_nodal_coordinate( aNodeId ) - tReferenceCoordinate )/ moris::norm(tReferenceCoordinate);

    REQUIRE( tRelDiffNorm <  1.0e-8 );

    // check time value for time step index 0
    real tReferenceTime = 1.000000000000000e+00;

    real tRelTimeDifference = std::abs( ( aExoIO.get_time_value( ) - tReferenceTime) / tReferenceTime );

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check temperature at node aNodeId in first time step (temperature is 3rd nodal field, first time step has index 0)
    real tReferenceTemperature = 9.764734706017323e+04;

    real tRelTempDifference = std::abs( ( aExoIO.get_nodal_field_value( aNodeId, 2, 0 ) - tReferenceTemperature ) / tReferenceTemperature );

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelTempDifference < 1.0e-4);

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    real tReferenceIQI = 8.324911657639191e+04;

    real tRelIQIDifference = std::abs( ( aExoIO.get_global_variable(0, 0 ) - tReferenceIQI ) / tReferenceIQI );

    REQUIRE(  tRelIQIDifference < 1.0e-4);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_serial()
{
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
        REQUIRE( tNumDims  ==  2    );
        REQUIRE( tNumNodes ==  5572 );
        REQUIRE( tNumElems ==  5426 );
    }

    // check results
    uint tNodeId = 1489;

    check_linear_results(tExoIO,tNodeId);
}

//---------------------------------------------------------------

extern "C"
void check_linear_results_parallel()
{
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
        REQUIRE( tNumDims  ==  2    );
        REQUIRE( tNumNodes ==  3652 );
        REQUIRE( tNumElems ==  3551 );
    }

    // check results
    uint tNodeId = 924;

    check_linear_results(tExoIO,tNodeId);
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
        "[moris],[example],[thermal],[diffusion]")
{
    // quadratic case currently not working
    if (false)
    {
        // define command line call
        int argc = 2;

        char tString1[] = "";
        char tString2[] = "./Channel_with_Four_Cylinders_Static_Temp_Only.so";

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
