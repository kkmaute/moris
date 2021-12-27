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
    Cell<uint> tReferenceNodeId  = {7805,2694,6977,4330};

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
    Cell<uint> tReferenceNumDims  = {3,3,3,3};
    Cell<uint> tReferenceNumNodes = {18292,6217,17824,5716};
    Cell<uint> tReferenceNumElems = {40254,13396,39168,12288};


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

    tReferenceCoordinate.push_back( { {+4.235886968816469e-01},{ +1.500000000000000e-01},{+2.358869688164689e-02} } );
    tReferenceCoordinate.push_back( { {+2.000000000000000e-01},{ +4.030968638358913e-01},{+0.000000000000000e+00} } );
    tReferenceCoordinate.push_back( { {+1.875000000000000e-01},{ +7.499999999999998e-02},{+3.875000000000000e-01} } );
    tReferenceCoordinate.push_back( { {+6.671327177817520e-02},{ +4.167132717781751e-01},{+4.999999999999999e-02} } );

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

    real tActualTime = tExoIO.get_time_value( );

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime(aTestCaseIndex)) / tReferenceTime(aTestCaseIndex) );

    MORIS_LOG_INFO("Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime(aTestCaseIndex),tActualTime,tRelDiffNorm*100.0);

    REQUIRE( tRelTimeDifference <  1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
     Cell<Matrix< DDRMat >> tReferenceDisplacement;

     tReferenceDisplacement.push_back( { {-2.695215210566065e+00},{-9.574313848576491e-01},{-1.509815634898354e-01} } );
     tReferenceDisplacement.push_back( { {-1.279576517421466e+00},{-2.572539354814230e+00},{-1.696472211536791e-05} } );
     tReferenceDisplacement.push_back( { {-1.265145383256113e+00},{-5.058657674628515e-01},{-2.567315492891605e+00} } );
     tReferenceDisplacement.push_back( { {-4.654107390524963e-01},{-2.880792909419323e+00},{-3.494429802007085e-01} } );

     Matrix< DDRMat > tActualDisplacement = {
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 2, 0 ) },
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 3, 0 ) },
             { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) } };

     real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement(aTestCaseIndex) ) / norm( tReferenceDisplacement(aTestCaseIndex) );

     MORIS_LOG_INFO("Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
             norm(tReferenceDisplacement(aTestCaseIndex)) ,norm(tActualDisplacement),tRelDispDifference*100.0);

    // FIXME: the displacement check is still failing for the "Immeresed" case
    REQUIRE(  tRelDispDifference < 1.0e-4);

    // check temperature at node aNodeId in first time step (temperature is 6th nodal field, first time step has index 0)
    Cell<real> tReferenceTemperature;
    tReferenceTemperature.push_back( 1.998065371606413e+02 );
    tReferenceTemperature.push_back( 1.996956411990300e+02 );
    tReferenceTemperature.push_back( 1.485698250771404e+02 );
    tReferenceTemperature.push_back( 9.988642380706278e+01 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 5, 0 );

    real tRelTempDifference = std::abs( ( tActualTemperature - tReferenceTemperature(aTestCaseIndex) ) / tReferenceTemperature(aTestCaseIndex) );

    MORIS_LOG_INFO("Check nodal temperature:    reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTemperature(aTestCaseIndex),tActualTemperature,tRelTempDifference*100.0);

    //FIXME: difference between parallel and serial run requires loose tolerance
    REQUIRE(  tRelTempDifference < 1.0e-4);

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Cell<real> tReferenceIQI;
    tReferenceIQI.push_back( 4.911266121905375e-01 );
    tReferenceIQI.push_back( 4.911266121905526e-01 );
    tReferenceIQI.push_back( 4.911031504660229e-01 );
    tReferenceIQI.push_back( 4.908693345404075e-01 );

    real tActualIQI = tExoIO.get_global_variable(0, 0 );

    real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI(aTestCaseIndex) ) / tReferenceIQI(aTestCaseIndex) );

    MORIS_LOG_INFO("Check temperature IQI:      reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceIQI(aTestCaseIndex),tActualIQI,tRelIQIDifference*100.0);

    REQUIRE(  tRelIQIDifference < 1.0e-4);
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
