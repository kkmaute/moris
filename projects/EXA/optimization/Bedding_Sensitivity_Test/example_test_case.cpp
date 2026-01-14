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
uint gLevelSetInterpolationOrder;
uint gFEMInterpolationOrder;
uint gLagrMeshInterpolationOrder;

// problem dimension: 2D or 3D
uint gDim;

// test case index
uint gTestCaseIndex;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//------------------------------------------------------------------------------------

void test_pause()
{
    // communicate process ids
    pid_t tPId = getpid();

    moris::Matrix< moris::DDSMat > tPIdVec;
    allgather_scalar( tPId, tPIdVec );

    if ( moris::par_rank() == 0 )
    {
        // print process Ids
        for ( int i = 0; i < moris::par_size(); ++i )
        {
            fprintf( stderr, "Process Rank %d ID: %d\n", i, tPIdVec( i ) );
        }

        std::string dummy;

        std::cout << "Press enter to continue . . .\n"
                  << std::flush;
        std::getline( std::cin, dummy );
    }
}

//---------------------------------------------------------------

extern "C" void
check_results(
        const std::string& aExoFileName,
        const std::string& aHdf5FileName,
        uint               aTestCaseIndex )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processors.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    //    // open and query exodus output file (set verbose to true to get basic mesh information)
    //    moris::mtk::Exodus_IO_Helper tExoIO(aExoFileName.c_str(),0,false,false);
    //
    //    // define reference node IDs
    //    Vector<uint> tReferenceNodeId  = {79,613,28,85};
    //
    //    // perturbation of denominator when building relative error
    //    real tDeltaEps = 1.0e-14;
    //
    //    if (gPrintReferenceValues)
    //    {
    //        std::cout << "Test case index: " << aTestCaseIndex << std::endl;
    //
    //        uint tNumDims  = tExoIO.get_number_of_dimensions();
    //        uint tNumNodes = tExoIO.get_number_of_nodes();
    //        uint tNumElems = tExoIO.get_number_of_elements();
    //
    //        std::cout << "Number of dimensions: " << tNumDims  << std::endl;
    //        std::cout << "Number of nodes     : " << tNumNodes << std::endl;
    //        std::cout << "Number of elements  : " << tNumElems << std::endl;
    //
    //        // coordinates of reference point
    //        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId(aTestCaseIndex) ), "Coordinates of reference point");
    //
    //        // time value for reference time step
    //        std::cout << "Time value: " << std::scientific << std::setprecision(15) << tExoIO.get_time_value() << std::endl;
    //
    //        // solution of reference point at reference time step
    //        std::cout << "Velocity in x-dir at reference point: " << std::scientific << std::setprecision(15) <<
    //                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 2, 0 ) << std::endl;
    //
    //        // solution of reference point at reference time step
    //        std::cout << "Velocity in y-dir at reference point: " << std::scientific << std::setprecision(15) <<
    //                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 3, 0 ) << std::endl;
    //
    //        // solution of reference point at reference time step
    //        std::cout << "Pressure at reference point: " << std::scientific << std::setprecision(15) <<
    //                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 ) << std::endl;
    //
    //        // solution of reference point at reference time step
    //        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision(15) <<
    //                tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 5, 0 ) << std::endl;
    //
    //        // value of IQI at reference time step
    //        hid_t tFileID = open_hdf5_file( aHdf5FileName );
    //        herr_t tStatus = 0;
    //
    //        Matrix<DDRMat> tConstraints;
    //        load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus);
    //
    //        std::cout << "IQI 0 value: " << std::scientific << std::setprecision(15) << tConstraints(0, 0 ) << std::endl;
    //        std::cout << "IQI 1 value: " << std::scientific << std::setprecision(15) << tConstraints(1, 0 ) << std::endl;
    //        std::cout << "IQI 2 value: " << std::scientific << std::setprecision(15) << tConstraints(2, 0 ) << std::endl;
    //        std::cout << "IQI 3 value: " << std::scientific << std::setprecision(15) << tConstraints(3, 0 ) << std::endl;
    //        std::cout << "IQI 4 value: " << std::scientific << std::setprecision(15) << tConstraints(4, 0 ) << std::endl;
    //        std::cout << "IQI 5 value: " << std::scientific << std::setprecision(15) << tConstraints(5, 0 ) << std::endl;
    //        std::cout << "IQI 6 value: " << std::scientific << std::setprecision(15) << tConstraints(6, 0 ) << std::endl;
    //        std::cout << "IQI 7 value: " << std::scientific << std::setprecision(15) << tConstraints(7, 0 ) << std::endl;
    //        std::cout << "IQI 8 value: " << std::scientific << std::setprecision(15) << tConstraints(8, 0 ) << std::endl;
    //        return;
    //    }
    //
    //    // define reference values for dimension, number of nodes and number of elements
    //    Vector<uint> tReferenceNumDims  = {  2,     2,    2,   2};
    //    Vector<uint> tReferenceNumNodes = { 654, 1043,  296, 413};
    //    Vector<uint> tReferenceNumElems = { 443,  443,  195, 195};
    //
    //    // check dimension, number of nodes and number of elements
    //    uint tNumDims  = tExoIO.get_number_of_dimensions();
    //    uint tNumNodes = tExoIO.get_number_of_nodes();
    //    uint tNumElems = tExoIO.get_number_of_elements();
    //
    //    MORIS_LOG_INFO("Check number of dimensions: reference %12d, actual %12d, percent  error %12.5e.",
    //            tReferenceNumDims(aTestCaseIndex),tNumDims,std::abs((tNumDims-tReferenceNumDims(aTestCaseIndex))/tReferenceNumDims(aTestCaseIndex)*100.0));
    //    MORIS_LOG_INFO("Check number of nodes:      reference %12d, actual %12d, percent  error %12.5e.",
    //            tReferenceNumNodes(aTestCaseIndex),tNumNodes,std::abs((tNumNodes-tReferenceNumNodes(aTestCaseIndex))/tReferenceNumNodes(aTestCaseIndex)*100.0));
    //    MORIS_LOG_INFO("Check number of elements:   reference %12d, actual %12d, percent  error %12.5e.",
    //            tReferenceNumElems(aTestCaseIndex),tNumElems,std::abs((tNumElems-tReferenceNumElems(aTestCaseIndex))/tReferenceNumElems(aTestCaseIndex)*100.0));
    //
    //    REQUIRE( tNumDims  ==  tReferenceNumDims(aTestCaseIndex)  );
    //    REQUIRE( tNumNodes ==  tReferenceNumNodes(aTestCaseIndex) );
    //    REQUIRE( tNumElems ==  tReferenceNumElems(aTestCaseIndex) );
    //
    //    // define reference coordinates for node aNodeId
    //    Vector<Matrix< DDRMat >> tReferenceCoordinate;
    //
    //    tReferenceCoordinate.push_back( { {1.041666666666667e+00},{5.892857142857143e-01} } );
    //    tReferenceCoordinate.push_back( { {1.041666666666667e+00},{5.892857142857142e-01} } );
    //    tReferenceCoordinate.push_back( { {1.041666666666667e+00},{5.892857142857143e-01} } );
    //    tReferenceCoordinate.push_back( { {1.041666666666667e+00},{5.892857142857143e-01} } );
    //
    //    // check nodal coordinates
    //    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId(aTestCaseIndex) );
    //
    //    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate(aTestCaseIndex) )/ (moris::norm(tReferenceCoordinate(aTestCaseIndex))+tDeltaEps);
    //
    //    MORIS_LOG_INFO("Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            tReferenceCoordinate(aTestCaseIndex)(0),tActualCoordinate(0),tRelDiffNorm*100.0);
    //    MORIS_LOG_INFO("Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            tReferenceCoordinate(aTestCaseIndex)(1),tActualCoordinate(1),tRelDiffNorm*100.0);
    //
    //    REQUIRE( tRelDiffNorm <  1.0e-8 );
    //
    //    // check time value for time step index 0
    //    Vector<real> tReferenceTime;
    //    tReferenceTime.push_back( 1.100000000000000e+01 );
    //    tReferenceTime.push_back( 1.100000000000000e+01 );
    //    tReferenceTime.push_back( 1.100000000000000e+01 );
    //    tReferenceTime.push_back( 1.100000000000000e+01 );
    //
    //    real tActualTime = tExoIO.get_time_value( );
    //
    //    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime(aTestCaseIndex)) / tReferenceTime(aTestCaseIndex) );
    //
    //    MORIS_LOG_INFO("Check time:                 reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            tReferenceTime(aTestCaseIndex),tActualTime,tRelTimeDifference*100.0);
    //
    //    REQUIRE( tRelTimeDifference <  1.0e-8 );
    //
    //    // check velocity at node aNodeId in first time step (velocities are 3rd and 4th  nodal field, first time step has index 0)
    //    Vector<Matrix<DDRMat>> tReferenceVelocity;
    //    tReferenceVelocity.push_back( { { 1.773958604669590e-01 } , {-5.709041615077051e-02} } );
    //    tReferenceVelocity.push_back( { { 1.082567779110537e-01 } , {-2.749907312681069e-02} } );
    //    tReferenceVelocity.push_back( { { 1.773958604669590e-01 } , {-5.709041615077050e-02} } );
    //    tReferenceVelocity.push_back( { { 1.082567779110537e-01 } , {-2.749907312681069e-02} } );
    //
    //    Matrix<DDRMat> tActualVelocity = {
    //            { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 2, 0 ) } ,
    //            { tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 3, 0 ) } };
    //
    //    real tRelVelocityDifference = norm( ( tActualVelocity - tReferenceVelocity(aTestCaseIndex) ) / norm ( tReferenceVelocity(aTestCaseIndex)+tDeltaEps ) );
    //
    //    MORIS_LOG_INFO("Check nodal velocities:     reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            norm ( tReferenceVelocity(aTestCaseIndex) ), norm( tActualVelocity ), tRelVelocityDifference*100.0);
    //
    //    REQUIRE(  tRelVelocityDifference < 1.0e-4);
    //
    //    // check pressure at node aNodeId in first time step (pressure is 5th nodal field, first time step has index 0)
    //    Vector<real> tReferencePressure;
    //    tReferencePressure.push_back( 8.794710548564350e-02 );
    //    tReferencePressure.push_back( 6.690311721734198e-02 );
    //    tReferencePressure.push_back( 8.794710548564350e-02 );
    //    tReferencePressure.push_back( 6.690311721734198e-02 );
    //
    //    real tActualPressure = tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 4, 0 );
    //
    //    real tRelPressureDifference = std::abs( ( tActualPressure - tReferencePressure(aTestCaseIndex) ) / ( tReferencePressure(aTestCaseIndex)+tDeltaEps ) );
    //
    //    MORIS_LOG_INFO("Check nodal Pressure:       reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            tReferencePressure(aTestCaseIndex),tActualPressure,tRelPressureDifference*100.0);
    //
    //    REQUIRE(  tRelPressureDifference < 1.0e-4);
    //
    //    // check temperature at node aNodeId in first time step (temperature is 6th nodal field, first time step has index 0)
    //    Vector<real> tReferenceTemperature;
    //    tReferenceTemperature.push_back( 2.920352898471725e-01);
    //    tReferenceTemperature.push_back( 4.357995793855248e-01 );
    //    tReferenceTemperature.push_back( 2.920352898471723e-01 );
    //    tReferenceTemperature.push_back( 4.357995793855248e-01 );
    //
    //    real tActualTemperature = tExoIO.get_nodal_field_value( tReferenceNodeId(aTestCaseIndex), 5, 0 );
    //
    //    real tRelTemperatureDifference = std::abs( ( tActualTemperature - tReferenceTemperature(aTestCaseIndex) ) / ( tReferenceTemperature(aTestCaseIndex)+tDeltaEps ) );
    //
    //    MORIS_LOG_INFO("Check nodal temperature:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //            tReferenceTemperature(aTestCaseIndex),tActualTemperature,tRelTemperatureDifference*100.0);
    //
    //    REQUIRE(  tRelTemperatureDifference < 1.0e-4);
    //
    //    // Sweep HDF5 file
    //    hid_t tFileID = open_hdf5_file( aHdf5FileName );
    //    herr_t tStatus = 0;
    //
    //    // Dedfine matrices for constraints and constraint gradients
    //    Matrix<DDRMat> tConstraints;
    //    Matrix<DDRMat> tConstraintsAnalytical;
    //    Matrix<DDRMat> tConstraintsFD;
    //
    //    // Check constraint values
    //    Vector<Matrix<DDRMat>> tReferenceConstraints;
    //
    //    tReferenceConstraints.push_back( {
    //        {-3.768956983866792e-06},
    //        { 6.902428439637079e-02},
    //        { 1.063063017173417e+00},
    //        { 8.236955476415919e-02},
    //        { 1.390581184223731e+00},
    //        {-3.638384242178661e-01},
    //        { 3.638374895164053e-01},
    //        { 1.507119843596119e-01},
    //        { 1.522513931060486e-01} } );
    //    tReferenceConstraints.push_back( {
    //        { 2.713219464337506e-04},
    //        { 7.511230990048032e-02},
    //        { 1.029062111905515e+00},
    //        { 3.290817990699941e-02},
    //        { 1.390581184223731e+00},
    //        {-2.327651687959512e-01},
    //        { 2.327643357675198e-01},
    //        { 1.497432653234213e-01},
    //        { 1.522513931060489e-01} } );
    //    tReferenceConstraints.push_back( {
    //        {-3.768956983866792e-06},
    //        { 6.902428439637079e-02},
    //        { 1.063063017173417e+00},
    //        { 8.236955476415919e-02},
    //        { 1.390581184223731e+00},
    //        {-3.638384242178661e-01},
    //        { 3.638374895164053e-01},
    //        { 1.507119843596119e-01},
    //        { 1.522513931060486e-01} } );
    //    tReferenceConstraints.push_back( {
    //        { 2.713219464337506e-04},
    //        { 7.511230990048032e-02},
    //        { 1.029062111905515e+00},
    //        { 3.290817990699941e-02},
    //        { 1.390581184223731e+00},
    //        {-2.327651687959512e-01},
    //        { 2.327643357675198e-01},
    //        { 1.497432653234213e-01},
    //        { 1.522513931060489e-01} } );
    //
    //    // Read constraints
    //    load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraints, tStatus);
    //
    //    for (uint tConIndex = 0; tConIndex < tConstraints.n_rows();++tConIndex)
    //    {
    //        real tRelConstraintDifference = std::abs( tConstraints(tConIndex) - tReferenceConstraints(aTestCaseIndex)(tConIndex) ) /
    //               ( std::abs( tReferenceConstraints(aTestCaseIndex)(tConIndex) ) +tDeltaEps );
    //
    //        MORIS_LOG_INFO("Check constraint %d:    reference %12.5e, actual %12.5e, percent  error %12.5e.",
    //                tReferenceConstraints(aTestCaseIndex)(tConIndex),tConstraints(tConIndex),tRelConstraintDifference*100.0);
    //
    //        REQUIRE(  tRelConstraintDifference < 1.0e-4);
    //    }

    // Sweep HDF5 file
    hid_t  tFileID = open_hdf5_file( aHdf5FileName );
    herr_t tStatus = 0;

    // Define matrices for constraints and constraint gradients
    Matrix< DDRMat > tConstraintsAnalytical;
    Matrix< DDRMat > tConstraintsFD;

    // Tolerance for adjoint vs. FD sensitivities
    moris::real tToleranceSensitivities = 1e-4;

    // perturbation of denominator when building relative error
    real tDeltaEps = 1.0e-12;

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

                if ( std::abs( tConstraintsFD( tConIndex, tADVIndex ) ) > 1e-6 )
                {
                    CHECK( tRelConstraintDifference < tToleranceSensitivities );
                }
                else
                {
                    CHECK( tConstraintsAnalytical( tConIndex, tADVIndex ) == Approx( tConstraintsFD( tConIndex, tADVIndex ) ).margin( 1e-6 ) );
                }
            }
        }
    }
}

//---------------------------------------------------------------

TEST_CASE( "Bedding_Sensitivity_Test",
        "[moris],[example],[optimization],[bedding]" )
{
    // remove files from previous test runs
    // FIXME: should be made independent of OS; note std::remove does not take wild cards
    if ( par_rank() == 0 )
    {
        MORIS_ERROR( std::system( "rm -f *exo*" ) == 0, "Bedding_Sensitivity_Test - removing *exo* files failed" );
        MORIS_ERROR( std::system( "rm -f *hdf5*" ) == 0, "Bedding_Sensitivity_Test - removing *hdf5* files failed" );
    }

    // for debugging in parallel
    // test_pause();

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Bedding_Sensitivity_Test.so";

    char* argv[ 2 ] = { tString1, tString2 };

    // Test cases 0 (serial) and 10 (parallel)

    // problem dimension: 2D or 3D
    gDim = 2;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Bedding Sensitivity Test - 2D - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    if ( par_size() == 1 )
    {
        // set test case index
        gTestCaseIndex = 0;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        barrier();

        // perform check for Test Case 0
        check_results( "Bedding_Sensitivity_Test_0.exo.e-s.0001", "Bedding_Sensitivity_Test_0_SEN.hdf5", gTestCaseIndex );
    }

    if ( par_size() == 4 )
    {
        // set test case index
        gTestCaseIndex = 10;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 2
        if ( par_rank() == 0 )
        {
            check_results( "Bedding_Sensitivity_Test_10.exo.e-s.0001", "Bedding_Sensitivity_Test_10_SEN.hdf5", gTestCaseIndex );
        }
    }

    // Test cases 1 (serial) and 11 (parallel)

    // problem dimension: 2D or 3D
    gDim = 3;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Bedding Sensitivity Test - 3D - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    if ( par_size() == 1 )
    {
        // set test case index
        gTestCaseIndex = 1;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        // perform check for Test Case 1
        check_results( "Bedding_Sensitivity_Test_1.exo.e-s.0001", "Bedding_Sensitivity_Test_1_SEN.hdf5", gTestCaseIndex );
    }

    if ( par_size() == 4 )
    {
        // set test case index
        gTestCaseIndex = 11;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        if ( par_rank() == 0 )
        {
            // perform check for Test Case 3
            check_results( "Bedding_Sensitivity_Test_11.exo.e-s.0001", "Bedding_Sensitivity_Test_11_SEN.hdf5", gTestCaseIndex );
        }
    }

#ifdef MORIS_HAVE_SLEPC
    // This cass tests the optimzation problemw with petsc solver
    if ( par_size() == 4 )
    {
        // set test case index
        gTestCaseIndex = 21;

        // call to performance manager main interface
        fn_WRK_Workflow_Main_Interface( argc, argv );

        if ( par_rank() == 0 )
        {
            // perform check for Test Case 3
            check_results( "Bedding_Sensitivity_Test_21.exo.e-s.0001", "Bedding_Sensitivity_Test_21_SEN.hdf5", gTestCaseIndex );
        }
    }
#endif
}
