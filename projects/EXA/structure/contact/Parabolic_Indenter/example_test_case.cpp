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
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "HDF5_Tools.hpp"

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

// text case index
uint gCaseIndex;

// Contact settings shared with the dynamically loaded nonlinear input.
bool gIndenterUseSeitz               = false;
bool gIndenterSensitivity            = false;
real gIndenterFrictionCoefficient    = -1.0;
real gIndenterTangentialDisplacement = 0.0;
uint gIndenterTimeSteps              = 1;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------

extern "C" void check_results()
{
    // check that run is serial; parallel not implemented yet
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" );

    // derive case index from linear/nonlinear string
    std::string tLinearNonlinear;

    switch ( gCaseIndex )
    {
        case 0:
            tLinearNonlinear = "Linear";
            break;
        case 1:
            tLinearNonlinear = "NonLinear";
            break;
        default:
            MORIS_ERROR( false, "check_results: gCaseIndex needs to be 0 or 1" );
    }

    std::string tExoFileName =
            "Parabolic_Indenter_" + tLinearNonlinear + "_Case_" + std::to_string( gCaseIndex ) + ".e-s.0000";

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", gCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( tExoFileName.c_str(), 0, false, false );

    // define reference node IDs
    // test case                          0    1
    Vector< uint > tReferenceNodeId = { 294, 104 };

    const uint tDispXIndex = 2;    // gCaseIndex = 0: IQIDispX, gCaseIndex = 1: IQIDispTopX
    const uint tDispYIndex = 3;    // gCaseIndex = 0: IQIDispY, gCaseIndex = 1: IQIDispTopY

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << gCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        // coordinates of reference point
        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId( gCaseIndex ) ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "Displacement at reference point: " << std::scientific << std::setprecision( 15 ) <<    //
                tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), tDispXIndex, 0 ) << "," <<     //
                tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), tDispYIndex, 0 ) << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    // test case                            0    1
    Vector< uint > tReferenceNumDims  = { 2, 2 };
    Vector< uint > tReferenceNumNodes = { 313, 313 };
    Vector< uint > tReferenceNumElems = { 192, 192 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims( gCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( gCaseIndex ) ) / tReferenceNumDims( gCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes( gCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( gCaseIndex ) ) / tReferenceNumNodes( gCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems( gCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( gCaseIndex ) ) / tReferenceNumElems( gCaseIndex ) * 100.0 ) );

    REQUIRE( tNumDims == tReferenceNumDims( gCaseIndex ) );
    REQUIRE( tNumNodes == tReferenceNumNodes( gCaseIndex ) );
    REQUIRE( tNumElems == tReferenceNumElems( gCaseIndex ) );

    // define reference coordinates for node aNodeId
    Vector< Matrix< DDRMat > > tReferenceCoordinate;

    tReferenceCoordinate.push_back( { { +5.000000000000001e-01 }, { +5.025000000000001e-01 } } );    // test case 0
    tReferenceCoordinate.push_back( { { +5.000000000000001e-01 }, { +5.034999999999999e-01 } } );    // test case 1

    // check nodal coordinates
    Matrix< DDRMat > tActualCoordinate = tExoIO.get_nodal_coordinate( tReferenceNodeId( gCaseIndex ) );

    real tRelDiffNorm = moris::norm( tActualCoordinate - tReferenceCoordinate( gCaseIndex ) ) / moris::norm( tReferenceCoordinate( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal x-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( gCaseIndex )( 0 ),
            tActualCoordinate( 0 ),
            std::abs( tActualCoordinate( 0 ) - tReferenceCoordinate( gCaseIndex )( 0 ) ) / moris::norm( tReferenceCoordinate( gCaseIndex ) ) * 100.0 );
    MORIS_LOG_INFO( "Check nodal y-coordinates:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceCoordinate( gCaseIndex )( 1 ),
            tActualCoordinate( 1 ),
            std::abs( tActualCoordinate( 1 ) - tReferenceCoordinate( gCaseIndex )( 1 ) ) / moris::norm( tReferenceCoordinate( gCaseIndex ) ) * 100.0 );

    REQUIRE( tRelDiffNorm < 1.0e-5 );

    // check time value for time step index 0
    Vector< real > tReferenceTime;
    tReferenceTime.push_back( 1.000000000000000e+01 );    // test case 0
    tReferenceTime.push_back( 1.000000000000000e+01 );    // test case 1

    real tActualTime = tExoIO.get_time_value();

    real tRelTimeDifference = std::abs( ( tActualTime - tReferenceTime( gCaseIndex ) ) / tReferenceTime( gCaseIndex ) );

    MORIS_LOG_INFO( "Check time:                 reference %12.5e, actual %12.5e, percent error %12.5e.",
            tReferenceTime( gCaseIndex ),
            tActualTime,
            tRelDiffNorm * 100.0 );

    REQUIRE( tRelTimeDifference < 1.0e-8 );

    // check displacements at node aNodeId in first time step (displacements are 3,4,5th nodal fields, first time step has index 0)
    Vector< Matrix< DDRMat > > tReferenceDisplacement;

    tReferenceDisplacement.push_back( { { -1.931735576379704e-05 }, { -2.458041207903132e-02 } } );    // test case 0
    tReferenceDisplacement.push_back( { { -1.159293612747665e-05 }, { -2.527716057948835e-02 } } );    // test case 1

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), tDispXIndex, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), tDispYIndex, 0 ) },
    };

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( gCaseIndex ) ) / norm( tReferenceDisplacement( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( gCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );
}

//---------------------------------------------------------------

TEST_CASE( "Parabolic_Indenter_Linear",
        "[moris],[example],[structure],[linear],[IndenterLegacy]" )
{
#ifdef MORIS_HAVE_ARBORX
    // check that run is serial; parallel not implemented yet
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Parabolic_Indenter_Linear.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    // set case index
    gCaseIndex = 0;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Parabolic_Indenter_Linear: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    check_results();

#else
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Parabolic_Indenter_Linear: Example skipped as Arborx not installed" );
    MORIS_LOG_INFO( " " );
#endif
}
//---------------------------------------------------------------

TEST_CASE( "Parabolic_Indenter_NonLinear",
        "[moris],[example],[structure],[nonlinear],[IndenterLegacy]" )
{
#ifdef MORIS_HAVE_ARBORX
    // check that run is serial; parallel not implemented yet
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" );

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Parabolic_Indenter_NonLinear.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // set interpolation order
    gInterpolationOrder = 1;

    // set case index
    gCaseIndex = 1;

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Parabolic_Indenter_NonLinear: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    check_results();

#else
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Parabolic_Indenter_NonLinear: Example skipped as Arborx not installed" );
    MORIS_LOG_INFO( " " );
#endif
}

#ifdef MORIS_HAVE_ARBORX
namespace
{
    void run_seitz( bool aFriction, bool aSensitivity )
    {
        REQUIRE( par_size() == 1 );
        gInterpolationOrder             = 1;
        gCaseIndex                      = 2 + ( aFriction ? 1 : 0 ) + ( aSensitivity ? 2 : 0 );
        gIndenterUseSeitz               = true;
        gIndenterSensitivity            = aSensitivity;
        gIndenterFrictionCoefficient    = aFriction ? 0.3 : -1.0;
        gIndenterTangentialDisplacement = 0.03;
        gIndenterTimeSteps              = 3;

        std::string tOutput = "Parabolic_Indenter_NonLinear_Case_" + std::to_string( gCaseIndex ) + ".exo";
        std::remove( tOutput.c_str() );
        std::remove( "Seitz_Indenter_Sensitivity.hdf5" );
        std::string         tLibrary = MORIS_INDENTER_NONLINEAR_LIBRARY;
        std::vector< char > tArgument( tLibrary.begin(), tLibrary.end() );
        tArgument.push_back( '\0' );
        char  tName[]      = "";
        char *tArguments[] = { tName, tArgument.data() };
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, tArguments ) == 0 );
    }

    void check_seitz_results()
    {
        const uint tLastStep = gIndenterTimeSteps - 1;
        // With Time_Offset=0, VIS uses a stable filename and physical time.
        const std::string tOutput = "Parabolic_Indenter_NonLinear_Case_" + std::to_string( gCaseIndex )
                                  + ".exo";
        mtk::Exodus_IO_Helper tExo( tOutput, tLastStep, false, false );
        REQUIRE( tExo.get_number_of_dimensions() == 2 );
        REQUIRE( tExo.get_number_of_nodes() == 313 );
        REQUIRE( tExo.get_number_of_elements() == 192 );
        const Matrix< DDRMat > tCoordinate = tExo.get_nodal_coordinate( 104 );
        CHECK( tCoordinate( 0 ) == Approx( 0.5 ).margin( 1e-12 ) );
        CHECK( tCoordinate( 1 ) == Approx( 0.5035 ).margin( 1e-12 ) );
        CHECK( tExo.get_time_value() == Approx( 10.0 ).margin( 1e-12 ) );

        // VIS fields 0/1 are node ID/index; displacement components follow.
        const real tUx        = tExo.get_nodal_field_value( 104, 2, tLastStep );
        const real tUy        = tExo.get_nodal_field_value( 104, 3, tLastStep );
        const real tYIntegral = tExo.get_global_variable( 0, tLastStep );
        const real tXIntegral = tExo.get_global_variable( 1, tLastStep );
        CHECK( std::isfinite( tUx ) );
        CHECK( std::isfinite( tUy ) );
        CHECK( std::isfinite( tYIntegral ) );
        CHECK( std::isfinite( tXIntegral ) );
        if ( gPrintReferenceValues )
        {
            std::cout << std::setprecision( 17 ) << "SEITZ_REFERENCE " << gCaseIndex << " "
                      << tUx << " " << tUy << " " << tYIntegral << " " << tXIntegral << '\n';
        }

        // Fresh serial three-step runs (2026-11-23): rows are frictionless
        // and mu=0.3; columns are UX, UY, integral(UY), integral(UX).
        // Sensitivity output follows the last small modulus perturbation, so
        // these field checks allow that change as well as solver roundoff.
        const real tReference[ 2 ][ 4 ] = {
            { 0.029768323798611617, -0.025033349010437068, -0.018893595366460098, 0.012356305254763569 },
            { 0.01701883707306719, -0.02520914297415712, -0.018904702175647648, 0.010933497998954718 }
        };
        const real tActual[] = { tUx, tUy, tYIntegral, tXIntegral };
        const uint tCase     = gIndenterFrictionCoefficient > 0.0 ? 1 : 0;
        for ( uint i = 0; i < 4; ++i )
        {
            CHECK( tActual[ i ] == Approx( tReference[ tCase ][ i ] ).epsilon( 1e-5 ).margin( 1e-8 ) );
        }
    }

    void check_seitz_sensitivities()
    {
        // The top modulus is the only ADV. Compare the final-state vertical
        // and horizontal displacement integrals against full workflow FDs.
        hid_t  tFile   = open_hdf5_file( "Seitz_Indenter_Sensitivity.hdf5" );
        herr_t tStatus = 0;
        for ( const std::string tQuantity : { "objective", "constraint" } )
        {
            Matrix< DDRMat > tAnalytical;
            load_matrix_from_hdf5_file( tFile, tQuantity + "_gradients eval_1-1 analytical", tAnalytical, tStatus );
            REQUIRE( tStatus == 0 );
            REQUIRE( tAnalytical.numel() == 1 );
            CHECK( std::isfinite( tAnalytical( 0 ) ) );
            CHECK( std::abs( tAnalytical( 0 ) ) > 1e-10 );
            for ( uint e = 1; e <= 2; ++e )
            {
                for ( const std::string tScheme : { "fd_forward", "fd_backward", "fd_central" } )
                {
                    Matrix< DDRMat > tFD;
                    load_matrix_from_hdf5_file( tFile, tQuantity + "_gradients eval_1-1 epsilon_" + std::to_string( e ) + "-2 " + tScheme, tFD, tStatus );
                    REQUIRE( tStatus == 0 );
                    REQUIRE( tFD.numel() == 1 );
                    INFO( tQuantity << ", epsilon " << e << ", " << tScheme );
                    CHECK( tAnalytical( 0 ) == Approx( tFD( 0 ) ).epsilon( 2e-4 ).margin( 1e-12 ) );
                    if ( gPrintReferenceValues )
                    {
                        std::cout << std::setprecision( 17 ) << "SEITZ_GRADIENT " << tQuantity << " "
                                  << tAnalytical( 0 ) << " " << tFD( 0 ) << '\n';
                    }
                }
            }
        }
        REQUIRE( close_hdf5_file( tFile ) == 0 );
    }
}    // namespace

// Nonlinear normal contact plus prescribed horizontal translation, with zero
// friction. Three increments exercise previous-state geometry and remapping.
TEST_CASE( "Parabolic indenter Seitz frictionless", "[SeitzFrictionless]" )
{
    run_seitz( false, false );
    check_seitz_results();
}

// Identical mesh, material and motion with Coulomb friction. The sequence
// carries previous displacement/traction through three physical load steps.
TEST_CASE( "Parabolic indenter Seitz friction", "[SeitzFriction]" )
{
    run_seitz( true, false );
    check_seitz_results();
}

// Check both final displacement integrals against two full-workflow FD step
// sizes for a material ADV, including frictionless contact in the adjoint.
TEST_CASE( "Parabolic indenter Seitz frictionless sensitivity", "[SeitzFrictionlessSensitivity]" )
{
    run_seitz( false, true );
    check_seitz_results();
    check_seitz_sensitivities();
}

// Repeat with friction history active, so the backward solve must include
// previous-state contact Jacobians rather than differentiating a frozen trace.
TEST_CASE( "Parabolic indenter Seitz friction sensitivity", "[SeitzFrictionSensitivity]" )
{
    run_seitz( true, true );
    check_seitz_results();
    check_seitz_sensitivities();
}
#endif
