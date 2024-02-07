//
// example specific interface to moris
//

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
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, true, true );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 4403, 4579, 2375, 2360 };

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
        std::cout << "Velocity at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << "," << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << "," <<
                // tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) <<
                std::endl;

        // solution of reference point at reference time step
        std::cout << "Pressure at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) << std::endl;

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 0, 0 ) << std::endl;

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2 };
    Vector< uint > tReferenceNumNodes = { 7127 };
    Vector< uint > tReferenceNumElems = { 11454 };

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

    tReferenceCoordinate.push_back( { { 0.0274990890175104 }, { 0.00591698987409472 } } );
    tReferenceCoordinate.push_back( { { 0.0316400639712811 }, { 0.00252554262988269 } } );
    tReferenceCoordinate.push_back( { { 0.0171121470630169 }, { 0.00168732192832977 } } );
    tReferenceCoordinate.push_back( { { 0.0171291530132294 }, { 0.000778483285102993 } } );

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
    // MORIS_LOG_INFO("Check nodal z-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
    // tReferenceCoordinate(aTestCaseIndex)(2),tActualCoordinate(2),tRelDiffNorm*100.0);

    REQUIRE( tRelDiffNorm < 1.0e-5 );

    // check time value for time step index 0
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 0.005 );
    tReferenceTime.push_back( 0.005 );
    tReferenceTime.push_back( 0.005 );
    tReferenceTime.push_back( 0.005 );
    // tReferenceTime.push_back( 1.000000000000000e+00 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceVelocity;

    tReferenceVelocity.push_back( { { -3.14130625281294e-07 }, { -2.69800723469803e-07 } } );
    tReferenceVelocity.push_back( { { -3.12490194316204e-07 }, { -5.38959617315861e-08 } } );
    tReferenceVelocity.push_back( { { -0.000377607866770003 }, { -1.55600063123971e-05 } } );
    tReferenceVelocity.push_back( { { -9.40252172704124e-06 }, { -3.36212068473784e-06 } } );

    Matrix< DDRMat > tActualVelocity = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) }
    };

    //{ tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) } };
    // std::cout << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << std::endl;
    // std::cout << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << std::endl;

    real tRelVelDifference = norm( tActualVelocity - tReferenceVelocity( aTestCaseIndex ) ) / norm( tReferenceVelocity( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal velocity:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceVelocity( aTestCaseIndex ) ),
            norm( tActualVelocity ),
            tRelVelDifference * 100.0 );

    REQUIRE( tRelVelDifference < 1.0e-5 );

    // check pressure at node aNodeId in first time step (temperature is 6th nodal field, first time step has index 0)
    Vector< real > tReferencePressure;
    tReferencePressure.push_back( 0.185522355835879 );
    tReferencePressure.push_back( 0.151972186085483 );
    tReferencePressure.push_back( 0.745731048974504 );
    tReferencePressure.push_back( 0.766439736246741 );

    real tActualPressure = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 );

    real tRelPresDifference = std::abs( ( tActualPressure - tReferencePressure( aTestCaseIndex ) ) / tReferencePressure( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal pressure:    reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferencePressure( aTestCaseIndex ),
            tActualPressure,
            tRelPresDifference * 100.0 );

    REQUIRE( tRelPresDifference < 1.0e-5 );

    // check IQI of first time step (only 1 IQI is defined, first time step has index 0)
    Vector< real > tReferenceIQI;
    tReferenceIQI.push_back( 0.185522355835879 );
    tReferenceIQI.push_back( 0.151972186085483 );
    tReferenceIQI.push_back( 0.745731048974504 );
    tReferenceIQI.push_back( 0.766439736246741 );

    // real tActualIQI = tExoIO.get_global_variable(0, 4 );
    // std::cout<<tActualIQI<<std::endl;

    // real tRelIQIDifference = std::abs( ( tActualIQI - tReferenceIQI(aTestCaseIndex) ) / tReferenceIQI(aTestCaseIndex) );

    // MORIS_LOG_INFO("Check Pressure IQI:      reference %12.5e, actual %12.5e, percent error %12.5e.",
    //     tReferenceIQI(aTestCaseIndex),tActualIQI,tRelIQIDifference*100.0);

    // REQUIRE(  tRelIQIDifference < 1.0e-5);
}

TEST_CASE( "Oscillator",
        "[moris],[example],[fluid],[BodyFitted]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Oscillator.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    // int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    // REQUIRE( tRet ==  0 );

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Oscillator: Interpolation order 1 - %i Processors.", par_size() );
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
            check_results( "Oscillator.exo.e-s.0000", 0 );
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
                check_results( "Oscillator.exo.e-s.0000", 1 );

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
