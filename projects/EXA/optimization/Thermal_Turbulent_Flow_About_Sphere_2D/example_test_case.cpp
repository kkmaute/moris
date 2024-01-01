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
    Vector< uint > tReferenceNodeId = { 350, 1266, 119, 408 };

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
        std::cout << "PhiD at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Velocity in x-dir at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Velocity in y-dir at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 5, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Pressure at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Viscosity at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 7, 0 ) << std::endl;

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 8, 0 ) << std::endl;

        // value of IQI at reference time step
        hid_t  tFileID = open_hdf5_file( aHdf5FileName );
        herr_t tStatus = 0;

        Matrix< DDRMat > tConstraints;
        load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus );

        std::cout << "IQI 0 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 0, 0 ) << std::endl;
        std::cout << "IQI 1 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 1, 0 ) << std::endl;
        std::cout << "IQI 2 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 2, 0 ) << std::endl;
        std::cout << "IQI 3 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 3, 0 ) << std::endl;
        std::cout << "IQI 4 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 4, 0 ) << std::endl;
        std::cout << "IQI 5 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 5, 0 ) << std::endl;
        std::cout << "IQI 6 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 6, 0 ) << std::endl;
        std::cout << "IQI 7 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 7, 0 ) << std::endl;
        std::cout << "IQI 8 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 8, 0 ) << std::endl;
        std::cout << "IQI 9 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 9, 0 ) << std::endl;
        std::cout << "IQI 10 value: " << std::scientific << std::setprecision( 15 ) << tConstraints( 10, 0 ) << std::endl;
        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2, 2, 2 };
    Vector< uint > tReferenceNumNodes = { 1440, 2785, 637, 1058 };
    Vector< uint > tReferenceNumElems = { 1307, 1307, 544, 544 };

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
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { 1.658018867924528e+00 }, { 7.403846153846154e-01 } } );
    tReferenceCoordinate.push_back( { { 1.658018867924528e+00 }, { 7.403846153846154e-01 } } );
    tReferenceCoordinate.push_back( { { 1.658018867924528e+00 }, { 7.403846153846154e-01 } } );
    tReferenceCoordinate.push_back( { { 1.658018867924528e+00 }, { 7.403846153846154e-01 } } );

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
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 1.100000000000000e+01 );
    tReferenceTime.push_back( 1.100000000000000e+01 );
    tReferenceTime.push_back( 1.100000000000000e+01 );
    tReferenceTime.push_back( 1.100000000000000e+01 );

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( aTestCaseIndex ) ) / tReferenceTime( aTestCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTime( aTestCaseIndex ),
            tActualTime,
            tRelTimeDifference * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check theta and phid at node aNodeId in first time step (first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceHMVars;
    tReferenceHMVars.push_back( { { 5.195264482367742e-01 }, { 2.959962414684608e-01 } } );
    tReferenceHMVars.push_back( { { 5.222488999432393e-01 }, { 2.925047833120029e-01 } } );
    tReferenceHMVars.push_back( { { 5.195264482367742e-01 }, { 2.959962414684608e-01 } } );
    tReferenceHMVars.push_back( { { 5.222488999432393e-01 }, { 2.925047833120029e-01 } } );

    Matrix< DDRMat > tActualHMVars = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 3, 0 ) }
    };

    real tRelHMVarsDifference = norm( ( tActualHMVars - tReferenceHMVars( aTestCaseIndex ) ) / norm( tReferenceHMVars( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal HM variables:     reference %12.5e, actual %12.5e, percent  error %12.5e.",
            norm( tReferenceHMVars( aTestCaseIndex ) ),
            norm( tActualHMVars ),
            tRelHMVarsDifference * 100.0 );

    REQUIRE( tRelHMVarsDifference < 1.0e-4 );

    // check velocity at node aNodeId in first time step (velocities are 3rd and 4th  nodal field, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceVelocity;
    tReferenceVelocity.push_back( { { 7.260479347944515e-01 }, { -7.372823057907386e-02 } } );
    tReferenceVelocity.push_back( { { 6.452713763324531e-01 }, { -8.756346273818047e-02 } } );
    tReferenceVelocity.push_back( { { 7.260479347944515e-01 }, { -7.372823057907386e-02 } } );
    tReferenceVelocity.push_back( { { 6.452713763324531e-01 }, { -8.756346273818047e-02 } } );

    Matrix< DDRMat > tActualVelocity = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 4, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 5, 0 ) }
    };

    real tRelVelocityDifference = norm( ( tActualVelocity - tReferenceVelocity( aTestCaseIndex ) ) / norm( tReferenceVelocity( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal velocities:     reference %12.5e, actual %12.5e, percent  error %12.5e.",
            norm( tReferenceVelocity( aTestCaseIndex ) ),
            norm( tActualVelocity ),
            tRelVelocityDifference * 100.0 );

    REQUIRE( tRelVelocityDifference < 1.0e-4 );

    // check pressure at node aNodeId in first time step (pressure is 5th nodal field, first time step has index 0)
    Vector< real > tReferencePressure;
    tReferencePressure.push_back( 1.683133156747930e-01 );
    tReferencePressure.push_back( 1.476207773463765e-01 );
    tReferencePressure.push_back( 1.683133156747930e-01 );
    tReferencePressure.push_back( 1.476207773463765e-01 );

    real tActualPressure = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 6, 0 );

    real tRelPressureDifference = std::abs( ( tActualPressure - tReferencePressure( aTestCaseIndex ) ) / ( tReferencePressure( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal Pressure:       reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferencePressure( aTestCaseIndex ),
            tActualPressure,
            tRelPressureDifference * 100.0 );

    REQUIRE( tRelPressureDifference < 1.0e-4 );

    // check turb. viscosity at node aNodeId in first time step (viscosity is 8th nodal field, first time step has index 0)
    Vector< real > tReferenceViscosity;
    tReferenceViscosity.push_back( 1.822749286708866e-03 );
    tReferenceViscosity.push_back( 1.702761040220779e-03 );
    tReferenceViscosity.push_back( 1.822749286708866e-03 );
    tReferenceViscosity.push_back( 1.702761040220779e-03 );

    real tActualViscosity = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 7, 0 );

    real tRelViscosityDifference = std::abs( ( tActualViscosity - tReferenceViscosity( aTestCaseIndex ) ) / ( tReferenceViscosity( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal viscosity:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceViscosity( aTestCaseIndex ),
            tActualViscosity,
            tRelViscosityDifference * 100.0 );

    REQUIRE( tRelViscosityDifference < 1.0e-4 );

    // check temperature at node aNodeId in first time step (temperature is 9th nodal field, first time step has index 0)
    Vector< real > tReferenceTemperature;
    tReferenceTemperature.push_back( 1.143523440655430e-01 );
    tReferenceTemperature.push_back( 1.271058787062102e-01 );
    tReferenceTemperature.push_back( 1.143523440655430e-01 );
    tReferenceTemperature.push_back( 1.271058787062102e-01 );

    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 8, 0 );

    real tRelTemperatureDifference = std::abs( ( tActualTemperature - tReferenceTemperature( aTestCaseIndex ) ) / ( tReferenceTemperature( aTestCaseIndex ) + tDeltaEps ) );

    MORIS_LOG_INFO( "Check nodal temperature:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
            tReferenceTemperature( aTestCaseIndex ),
            tActualTemperature,
            tRelTemperatureDifference * 100.0 );

    REQUIRE( tRelTemperatureDifference < 1.0e-4 );

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( aHdf5FileName );
    herr_t tStatus = 0;

    // Define matrices for constraints and constraint gradients
    Matrix< DDRMat > tConstraints;
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tConstraintsFD;

    // Tolerance for adjoint vs. FD sensitivities
    moris::real tToleranceSensivities = 0.001;

    // Check constraint values
    Vector< Matrix< DDRMat > > tReferenceConstraints;

    tReferenceConstraints.push_back( { //
            { 2.298131639762928e-01 },
            { 1.107447112250336e+00 },
            { 1.324837464112712e-01 },
            { 1.407259385555996e+00 },
            { -4.677853462592048e-01 },
            { 4.677853462592049e-01 },
            { 1.545543281610846e-01 },
            { 1.565261379254385e-01 },
            { 5.241036081238704e-01 },
            { -7.908056694053146e-02 } } );

    tReferenceConstraints.push_back( { //
            { 2.218451925583727e-01 },
            { 1.091344806154332e+00 },
            { 1.057239328238967e-01 },
            { 1.407259385555996e+00 },
            { -4.184211403842039e-01 },
            { 4.184211403842040e-01 },
            { 1.542375136588137e-01 },
            { 1.565261379254387e-01 },
            { 4.627550584962253e-01 },
            { -5.629465456553837e-02 } } );

    tReferenceConstraints.push_back( { //
            { 2.298131639762928e-01 },
            { 1.107447112250336e+00 },
            { 1.324837464112712e-01 },
            { 1.407259385555996e+00 },
            { -4.677853462592048e-01 },
            { 4.677853462592049e-01 },
            { 1.545543281610846e-01 },
            { 1.565261379254385e-01 },
            { 5.241036081238704e-01 },
            { -7.908056694053146e-02 } } );

    tReferenceConstraints.push_back( { //
            { 2.218451925583727e-01 },
            { 1.091344806154332e+00 },
            { 1.057239328238967e-01 },
            { 1.407259385555996e+00 },
            { -4.184211403842039e-01 },
            { 4.184211403842040e-01 },
            { 1.542375136588137e-01 },
            { 1.565261379254387e-01 },
            { 4.627550584962253e-01 },
            { -5.629465456553837e-02 } } );

    // Read constraints
    load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus );

    for ( uint tConIndex = 0; tConIndex < tConstraints.n_rows(); ++tConIndex )
    {
        real tRelConstraintDifference = std::abs( tConstraints( tConIndex ) - tReferenceConstraints( aTestCaseIndex )( tConIndex ) ) /    //
                                        ( std::abs( tReferenceConstraints( aTestCaseIndex )( tConIndex ) ) + tDeltaEps );

        MORIS_LOG_INFO( "Check constraint %d:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
                tConIndex,
                tReferenceConstraints( aTestCaseIndex )( tConIndex ),
                tConstraints( tConIndex ),
                100.0 * tRelConstraintDifference );

        REQUIRE( tRelConstraintDifference < 1.0e-4 );
    }

    // Read analytical sensitivities
    load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintsAnalytical, tStatus );

    // Read FD sensitivities and compare
    Vector< std::string > tFDTypes = { "fd_central" };
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

                if ( std::abs( tConstraintsFD( tADVIndex ) ) > 1e-6 )
                {
                    CHECK( tRelConstraintDifference < tToleranceSensivities );
                }
                else
                {
                    CHECK( tConstraintsAnalytical( tADVIndex ) == Approx( tConstraintsFD( tADVIndex ) ).margin( 1e-8 ) );
                }
            }
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Thermal_Turbulent_Flow_About_Sphere_2D",
        "[moris],[example],[optimization],[advection],[turbulence]" )
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./ChannelTurbulent2D.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Thermal Turbulent Flow about Sphere - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    if ( par_size() == 1 )
    {
        // set test case index
        gTestCaseIndex = 0;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 0
        check_results( "ChannelTurbulent2D_0.exo.e-s.0001", "ChannelTurbulent2D_0_SEN.hdf5", gTestCaseIndex );
    }

    if ( par_size() == 3 )
    {
        // set test case index
        gTestCaseIndex = 2;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 2
        if ( par_rank() == 0 )
        {
            check_results( "ChannelTurbulent2D_2.exo.e-s.0001", "ChannelTurbulent2D_2_SEN.hdf5", gTestCaseIndex );
        }
    }

    // set interpolation order
    gInterpolationOrder = 2;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Thermal Turbulent Flow about Sphere - 2D: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    if ( par_size() == 1 )
    {
        // set test case index
        gTestCaseIndex = 1;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 1
        check_results( "ChannelTurbulent2D_1.exo.e-s.0001", "ChannelTurbulent2D_1_SEN.hdf5", gTestCaseIndex );
    }

    if ( par_size() == 3 )
    {
        // set test case index
        gTestCaseIndex = 3;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        if ( par_rank() == 0 )
        {
            // perform check for Test Case 3
            check_results( "ChannelTurbulent2D_3.exo.e-s.0001", "ChannelTurbulent2D_3_SEN.hdf5", gTestCaseIndex );
        }
    }
}
