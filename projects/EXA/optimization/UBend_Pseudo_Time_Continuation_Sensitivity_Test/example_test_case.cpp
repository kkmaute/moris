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

// test case index
uint gTestCaseIndex;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

extern "C" void
check_results(
        std::string aExoFileName,
        std::string aHdf5FileName,
        uint        aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processors.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    Cell< uint > tReferenceNodeId = { 1991, 735 };

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
        std::cout << "Velocity in x-dir at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 5, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Velocity in y-dir at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Pressure at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 7, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Viscosity at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 8, 0 ) << std::endl;

        // value of IQI at reference time step
        hid_t  tFileID = open_hdf5_file( aHdf5FileName );
        herr_t tStatus = 0;

        Matrix< DDRMat > tConstraints;
        load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus );

        std::cout << "IQI 0 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 0 ) << std::endl;
        std::cout << "IQI 1 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 1 ) << std::endl;
        std::cout << "IQI 2 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 2 ) << std::endl;
        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Cell< uint > tReferenceNumDims  = { 2, 2 };
    Cell< uint > tReferenceNumNodes = { 2916, 848 };
    Cell< uint > tReferenceNumElems = { 2607, 772 };

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

    tReferenceCoordinate.push_back( { { 6.991440000000001e+00 }, { 3.486625000000000e+00 } } );
    tReferenceCoordinate.push_back( { { 6.991440000000001e+00 }, { 3.486625000000000e+00 } } );

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

    REQUIRE( tRelDiffNorm < 1.0e-8 );

    // check time value for time step index 0
    Cell< real > tReferenceTime;
    tReferenceTime.push_back( 1.100000000000000e+01 );
    tReferenceTime.push_back( 1.100000000000000e+01 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelTimeDifference * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check velocity at node aNodeId in first time step (velocities are 3rd and 4th  nodal field, first time step has index 0)
    Cell< Matrix< DDRMat > > tReferenceVelocity;
    tReferenceVelocity.push_back( { { -2.621985520118972e-01 }, { -1.223821117360914e-01 } } );
    tReferenceVelocity.push_back( { { -2.621994928927610e-01 }, { -1.223827050341624e-01 } } );

    Matrix< DDRMat > tActualVelocity = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 5, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 ) }
    };

    real tRelVelocityDifference = norm( ( tActualVelocity - tReferenceVelocity( aTestCaseIndex ) ) / norm( tReferenceVelocity( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal velocities:     reference %12.5e, actual %12.5e, percent  error %12.5e.",
            norm( tReferenceVelocity( aTestCaseIndex ) ),
            norm( tActualVelocity ),
            tRelVelocityDifference * 100.0 );

    REQUIRE( tRelVelocityDifference < 1.0e-4 );

    // check pressure at node aNodeId in first time step (pressure is 5th nodal field, first time step has index 0)
    Cell< real > tReferencePressure;
    tReferencePressure.push_back( 6.600450821884006e-01 );
    tReferencePressure.push_back( 6.600454151170850e-01 );

    real tActualPressure = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 7, 0 );

    real tRelPressureDifference = std::abs( ( tActualPressure - tReferencePressure( aTestCaseIndex ) ) / ( tReferencePressure( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Pressure:       reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferencePressure( aTestCaseIndex ),
            tActualPressure,
            tRelPressureDifference * 100.0 );

    REQUIRE( tRelPressureDifference < 1.0e-4 );

    // check Viscosity at node aNodeId in first time step (Viscosity is 6th nodal field, first time step has index 0)
    Cell< real > tReferenceViscosity;
    tReferenceViscosity.push_back( 1.967343367455367e-02 );
    tReferenceViscosity.push_back( 1.967348194097030e-02 );

    real tActualViscosity = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 8, 0 );

    real tRelViscosityDifference = std::abs( ( tActualViscosity - tReferenceViscosity( aTestCaseIndex ) ) / ( tReferenceViscosity( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Viscosity:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceViscosity( aTestCaseIndex ),
            tActualViscosity,
            tRelViscosityDifference * 100.0 );

    REQUIRE( tRelViscosityDifference < 1.0e-4 );

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( aHdf5FileName );
    herr_t tStatus = 0;

    // Dedfine matrices for constraints and constraint gradients
    Matrix< DDRMat > tConstraints;
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tConstraintsFD;

    // Tolerance for adjoint vs. FD sensitivities
    moris::real tToleranceSensivities = 0.001;

    // Check constraint values
    Cell< Matrix< DDRMat > > tReferenceConstraints;

    tReferenceConstraints.push_back( { { 1.158612227050382e+01 },
            { 2.120790687504068e+01 },
            { 1.023252485864117e+02 } } );
    tReferenceConstraints.push_back( { { 1.158612227050381e+01 },
            { 2.120790687504071e+01 },
            { 1.023252485864111e+02 } } );

    // Read constraints
    load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus );

    for ( uint tConIndex = 0; tConIndex < tConstraints.n_rows(); ++tConIndex )
    {
        real tRelConstraintDifference = std::abs( tConstraints( tConIndex ) - tReferenceConstraints( aTestCaseIndex )( tConIndex ) ) / ( std::abs( tReferenceConstraints( aTestCaseIndex )( tConIndex ) ) + tDeltaEps );

        MORIS_LOG_INFO( "Check constraint %d:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
                tConIndex,
                tReferenceConstraints( aTestCaseIndex )( tConIndex ),
                tConstraints( tConIndex ),
                tRelConstraintDifference * 100.0 );

        REQUIRE( tRelConstraintDifference < 1.0e-4 );
    }

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );

    // Read FD sensitivities and compare
    Cell< std::string > tFDTypes = { "fd_forward", "fd_backward", "fd_central" };
    for ( uint tFDIndex = 0; tFDIndex < tFDTypes.size(); tFDIndex++ )
    {
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 " + tFDTypes( tFDIndex ), tConstraintsFD, tStatus );

        REQUIRE( tConstraintsAnalytical.numel() == tConstraintsFD.numel() );

        for ( uint tADVIndex = 0; tADVIndex < tConstraintsAnalytical.n_cols(); tADVIndex++ )
        {
            for ( uint tConIndex = 0; tConIndex < tConstraintsAnalytical.n_rows(); ++tConIndex )
            {
                real tRelConstraintDifference =
                        std::abs( tConstraintsAnalytical( tConIndex, tADVIndex ) - tConstraintsFD( tConIndex, tADVIndex ) ) /    //
                        ( std::abs( tConstraintsFD( tConIndex, tADVIndex ) ) + tDeltaEps );

                MORIS_LOG_INFO( "Check derivative of constraint %d  wrt. ADV(%i):  analytical  %12.5e, finite difference (%s) %12.5e, percent error %12.5e.",
                        tConIndex,
                        tADVIndex,
                        tConstraintsAnalytical( tConIndex, tADVIndex ),
                        tFDTypes( tFDIndex ).c_str(),
                        tConstraintsFD( tConIndex, tADVIndex ),
                        100.0 * tRelConstraintDifference );

                CHECK( tRelConstraintDifference < tToleranceSensivities );
            }
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "UBend_Pseudo_Time_Continuation_Sensitivity_Test",
        "[moris],[example],[optimization],[turbulence]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./uBend_2D_Optimization.so";

    char* argv[ 2 ] = { tString1, tString2 };

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "UBend_Pseudo_Time_Continuation_Sensitivity_Test:%i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    if ( par_size() == 1 )
    {
        // set test case index
        gTestCaseIndex = 0;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 0
        check_results( "uBend_2D_Optimization.exo.e-s.0001", "uBend_2D_Optimization.hdf5", gTestCaseIndex );
    }

    if ( par_size() == 4 )
    {
        // set test case index
        gTestCaseIndex = 1;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 2
        if ( par_rank() == 0 )
        {
            check_results( "uBend_2D_Optimization.exo.e-s.0001", "uBend_2D_Optimization.hdf5", gTestCaseIndex );
        }
    }
}