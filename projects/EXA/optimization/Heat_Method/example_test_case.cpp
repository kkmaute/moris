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

#include "HDF5_Tools.hpp"

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
        std::string aHdf5FileName,
        uint        aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Cell< uint > tReferenceNodeId = { 4, 13, 4, 13 };

    // perturbation of denominator when building relative error
    real tDeltaEps = 1.0e-14;

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
        std::cout << "Theta at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Phid at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Phi-design (Pdsg) at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) << std::endl;
        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell< uint > tReferenceNumDims  = { 2, 3, 2, 3 };
    Cell< uint > tReferenceNumNodes = { 82, 414, 130, 1088 };
    Cell< uint > tReferenceNumElems = { 60, 632, 60, 920 };

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
    Cell< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { 0.0 }, { 0.0 } } );
    tReferenceCoordinate.push_back( { { 0.0 }, { 0.0 }, { 0.0 } } );
    tReferenceCoordinate.push_back( { { 0.0 }, { 0.0 } } );
    tReferenceCoordinate.push_back( { { 0.0 }, { 0.0 }, { 0.0 } } );

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( aTestCaseIndex ) ) / ( moris::norm( tReferenceCoordinate( aTestCaseIndex ) ) + tDeltaEps );

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
    Cell< real > tReferenceTime;
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

    // check Theta at node aNodeId in first time step (Theta is 3rd nodal field, first time step has index 0)
    Cell< real > tReferenceTheta;
    tReferenceTheta.push_back( 5.917299018173376e-01 );
    tReferenceTheta.push_back( 7.194331939566655e-01 );
    tReferenceTheta.push_back( 6.349840630136200e-01 );
    tReferenceTheta.push_back( 7.508299324748982e-01 );

    real tActualTheta = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 );

    real tRelThetaDifference = std::abs( ( tActualTheta - tReferenceTheta( aTestCaseIndex ) ) / ( tReferenceTheta( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Theta:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTheta( aTestCaseIndex ),
            tActualTheta,
            tRelThetaDifference * 100.0 );

    REQUIRE( tRelThetaDifference < 1.0e-4 );

    // check Phid at node aNodeId in first time step (Phid is 4th nodal field, first time step has index 0)
    Cell< real > tReferencePhid;
    tReferencePhid.push_back( 3.872167974340174e-01 );
    tReferencePhid.push_back( 3.587919624293675e-01 );
    tReferencePhid.push_back( 2.740283249190941e-01 );
    tReferencePhid.push_back( 2.459025775192394e-01 );

    real tActualPhid = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 );

    real tRelPhidDifference = std::abs( ( tActualPhid - tReferencePhid( aTestCaseIndex ) ) / ( tReferencePhid( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Phid:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferencePhid( aTestCaseIndex ),
            tActualPhid,
            tRelPhidDifference * 100.0 );

    REQUIRE( tRelPhidDifference < 1.0e-4 );

    // check Pdsg at node aNodeId in first time step (Pdsg is 4th nodal field, first time step has index 0)
    Cell< real > tReferencePdsg;
    tReferencePdsg.push_back( 4.710000000000003e-01 );
    tReferencePdsg.push_back( 4.709999999999997e-01 );
    tReferencePdsg.push_back( 2.862624741953033e-01 );
    tReferencePdsg.push_back( 2.506857164624124e-01 );

    real tActualPdsg = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 );

    real tRelPdsgDifference = std::abs( ( tActualPdsg - tReferencePdsg( aTestCaseIndex ) ) / ( tReferencePdsg( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Pdsg:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferencePdsg( aTestCaseIndex ),
            tActualPdsg,
            tRelPdsgDifference * 100.0 );

    REQUIRE( tRelPdsgDifference < 1.0e-4 );

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( aHdf5FileName );
    herr_t tStatus = 0;

    // Declare sensitivity matrices for comparison
    Matrix< DDRMat > tObjectiveAnalytical;
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tObjectiveFD;
    Matrix< DDRMat > tConstraintsFD;

    // Tolerance for adjoint vs. FD sensitivities
    moris::real tToleranceSensties = 0.001;

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytical, tStatus );
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );

    // Read FD sensitivities and compare
    Cell< std::string > tFDTypes = { "fd_forward", "fd_backward", "fd_central" };
    for ( uint tFDIndex = 0; tFDIndex < tFDTypes.size(); tFDIndex++ )
    {
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tObjectiveFD, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tConstraintsFD, tStatus );

        REQUIRE( tObjectiveAnalytical.numel() == tObjectiveFD.numel() );
        REQUIRE( tConstraintsAnalytical.numel() == tConstraintsFD.numel() );

        for ( uint tADVIndex = 0; tADVIndex < tObjectiveAnalytical.length(); tADVIndex++ )
        {
            MORIS_LOG_INFO( "Check derivative of objective  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                    tADVIndex,
                    tObjectiveAnalytical( tADVIndex ),
                    tFDTypes( tFDIndex ).c_str(),
                    tObjectiveFD( tADVIndex ),
                    100 * std::abs( ( tObjectiveAnalytical( tADVIndex ) - tObjectiveFD( tADVIndex ) ) / ( tObjectiveFD( tADVIndex ) + tDeltaEps ) ) );

            CHECK( std::abs( ( tObjectiveAnalytical( tADVIndex ) - tObjectiveFD( tADVIndex ) ) / ( tObjectiveFD( tADVIndex ) + tDeltaEps ) ) < tToleranceSensties );

            for ( uint tConIndex = 0; tConIndex < tConstraintsAnalytical.n_rows(); ++tConIndex )
            {
                MORIS_LOG_INFO( "Check derivative of constraint %d  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                        tConIndex,
                        tADVIndex,
                        tConstraintsAnalytical( tConIndex, tADVIndex ),
                        tFDTypes( tFDIndex ).c_str(),
                        tConstraintsFD( tConIndex, tADVIndex ),
                        100 * std::abs( ( tConstraintsAnalytical( tConIndex, tADVIndex ) - tConstraintsFD( tConIndex, tADVIndex ) ) / ( tConstraintsFD( tConIndex, tADVIndex ) + tDeltaEps ) ) );

                CHECK( std::abs( ( tConstraintsAnalytical( tConIndex, tADVIndex ) - tConstraintsFD( tConIndex, tADVIndex ) ) / ( tConstraintsFD( tConIndex, tADVIndex ) + tDeltaEps ) ) < tToleranceSensties );
            }
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "HeatMethod_Linear",
        "[moris],[example],[structure],[linear]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./HeatMethod.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Heat Method - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 0;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 0
    check_results( "HeatMethod_0.exo", "SEN_HeatMethod_0.hdf5", gTestCaseIndex );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Heat Method - 3D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 3D
    gDim = 3;

    // set test case index
    gTestCaseIndex = 1;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 1
    check_results( "HeatMethod_1.exo", "SEN_HeatMethod_1.hdf5", gTestCaseIndex );
}

//---------------------------------------------------------------

TEST_CASE( "HeatMethod_Quadratic",
        "[moris],[example],[structure],[quadratic]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./HeatMethod.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Heat Method - 2D: Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 2D
    gDim = 2;

    // set test case index
    gTestCaseIndex = 2;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 2
    check_results( "HeatMethod_2.exo", "SEN_HeatMethod_2.hdf5", gTestCaseIndex );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Heat Method - 3D: Interpolation order 2 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // set dimension: 3D
    gDim = 3;

    // set test case index
    gTestCaseIndex = 3;

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // perform check for Test Case 3
    check_results( "HeatMethod_3.exo", "SEN_HeatMethod_3.hdf5", gTestCaseIndex );
}
