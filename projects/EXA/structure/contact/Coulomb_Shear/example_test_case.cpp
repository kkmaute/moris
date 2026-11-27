/*
 * Copyright (c) 2026 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 */
#include "catch.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "moris_typedefs.hpp"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

moris::uint gShearElements        = 8;
moris::uint gShearSubsteps        = 1;
moris::real gShearFriction        = 0.3;
moris::real gShearYoungsModulus   = 1000.0;
moris::real gShearCompression     = 1.0e-4;
moris::real gShearPenalty         = 500.0;
bool        gPrintReferenceValues = false;

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

namespace
{
    using namespace moris;

    // Analytic endpoints in units of Dcrit and mu*p. These are the closed-form
    // two-layer series compliance and Coulomb slider solution, not solver output.
    const real     tDrive[] = { 0, 0.25, 0.5, 1, 1.5, 2, 2, 1.5, 1, 0.5, 0, -0.5, -2, -2, -1.5, -1, -0.5, 0, 0.5, 2 };
    const real     tShear[] = { 0, 0.25, 0.5, 1, 1, 1, 1, 0.5, 0, -0.5, -1, -1, -1, -1, -0.5, 0, 0.5, 1, 1, 1 };
    const real     tSlip[]  = { 0, 0, 0, 0, 0.5, 1, 1, 1, 1, 1, 1, 0.5, -1, -1, -1, -1, -1, -1, -0.5, 1 };
    constexpr uint tStages  = 20;

    uint reference_node( mtk::Exodus_IO_Helper& aOutput, real aY )
    {
        for ( uint i = 0; i < static_cast< uint >( aOutput.get_number_of_nodes() ); ++i )
        {
            const auto tId = aOutput.get_node_id( i );
            const auto tX  = aOutput.get_nodal_coordinate( tId );
            if ( std::abs( tX( 0 ) - 0.5 ) < 1.0e-10 && std::abs( tX( 1 ) - aY ) < 1.0e-10 )
            {
                return tId;
            }
        }
        FAIL( "Missing central reference node" );
        return 0;
    }

    void check_results()
    {
        // Serial tests: this rank owns both interior reference nodes. Keep any
        // collective operations outside the reader/checker lifetime.
        REQUIRE( par_size() == 1 );
        const real tG          = 0.5 * gShearYoungsModulus;
        const real tLambda     = 1.0 - gShearCompression / 0.98;
        const real tPressure   = tG * ( 1.0 / tLambda - tLambda );
        const real tLimit      = gShearFriction * tPressure;
        const real tDcrit      = tLimit * 0.98 / tG;
        uint       tTopNode    = 0;
        uint       tBottomNode = 0;
        {
            mtk::Exodus_IO_Helper tOutput( "Coulomb_Shear.exo", 0, false, false );
            REQUIRE( tOutput.get_number_of_dimensions() == 2 );
            // Serial Q1 topology references obtained on 2026-11-26. The
            // physical force/slip references below are analytical, not fitted.
            REQUIRE( tOutput.get_number_of_nodes() == ( gShearElements == 16 ? 470 : 174 ) );
            REQUIRE( tOutput.get_number_of_elements() == ( gShearElements == 16 ? 416 : 144 ) );
            tTopNode    = reference_node( tOutput, 0.75 );
            tBottomNode = reference_node( tOutput, 0.25 );
            CHECK( tTopNode == ( gShearElements == 16 ? 294u : 100u ) );
            CHECK( tBottomNode == ( gShearElements == 16 ? 76u : 22u ) );
            if ( gPrintReferenceValues )
            {
                std::cout << "COULOMB_MESH " << gShearElements << " " << tOutput.get_number_of_nodes()
                          << " " << tOutput.get_number_of_elements() << " " << tTopNode << " " << tBottomNode << std::endl;
            }
        }

        std::ofstream tHistory( "Coulomb_Shear_history.csv" );
        tHistory << "time,drive,shear,normal,slip,reference_shear,reference_slip\n"
                 << std::setprecision( 17 );
        real tLastSlip = 0.0;
        for ( uint i = 0; i < tStages; ++i )
        {
            const uint tStep = ( i + 1 ) * gShearSubsteps - 1;
            // Initialize both nodal and global caches at this step; older helper
            // versions share a cache index. Close each reader before reopening.
            mtk::Exodus_IO_Helper tOutput( "Coulomb_Shear.exo", tStep, false, false );
            CHECK( tOutput.get_time_value() == Approx( static_cast< real >( ( i + 1 ) * gShearSubsteps ) ).margin( 1e-12 ) );
            const real tShearForce  = tOutput.get_global_variable( 0, tStep );
            const real tNormalForce = -tOutput.get_global_variable( 1, tStep );
            // Fields 0/1 are VIS node metadata; each phase has its own UX/UY.
            const real tTopUx        = tOutput.get_nodal_field_value( tTopNode, 2, tStep );
            const real tBottomUx     = tOutput.get_nodal_field_value( tBottomNode, 4, tStep );
            const real tMeasuredSlip = 1.96 * ( tTopUx - tBottomUx ) - 0.96 * tDrive[ i ] * tDcrit;
            CAPTURE( i, tShearForce, tNormalForce, tMeasuredSlip );
            REQUIRE( std::isfinite( tShearForce ) );
            REQUIRE( std::isfinite( tNormalForce ) );
            REQUIRE( std::isfinite( tMeasuredSlip ) );
            if ( gPrintReferenceValues )
            {
                std::cout << std::setprecision( 15 ) << "COULOMB_STAGE " << i << " " << tShearForce
                          << " " << tNormalForce << " " << tMeasuredSlip << " ux " << tTopUx << " " << tBottomUx << std::endl;
            }
            CHECK( tNormalForce == Approx( tPressure ).epsilon( 0.005 ).margin( 1.0e-8 ) );
            CHECK( tShearForce == Approx( tShear[ i ] * tLimit ).margin( 0.02 * tLimit ) );
            CHECK( tMeasuredSlip == Approx( tSlip[ i ] * tDcrit ).margin( 0.01 * tDcrit ) );
            CHECK( tTopUx == Approx( tDrive[ i ] * tDcrit - 0.25 * tShear[ i ] * tLimit / tG ).margin( 0.005 * tDcrit ) );
            CHECK( tBottomUx == Approx( 0.25 * tShear[ i ] * tLimit / tG ).margin( 0.005 * tDcrit ) );
            const real tTopUy    = tOutput.get_nodal_field_value( tTopNode, 3, tStep );
            const real tBottomUy = tOutput.get_nodal_field_value( tBottomNode, 5, tStep );
            const real tGap      = 0.02 + 1.96 * ( tTopUy - tBottomUy ) + 0.96 * ( 0.02 + gShearCompression );
            CHECK( std::abs( tGap ) < 1.0e-7 );
            CHECK( std::abs( tShearForce ) <= gShearFriction * tNormalForce + 0.005 * tLimit );
            // The sign convention reports the resisting shear magnitude with
            // the imposed drive; shear times relative slip is dissipated work.
            CHECK( tShearForce * ( tMeasuredSlip - tLastSlip ) >= -0.005 * tLimit * tDcrit );
            tHistory << ( i + 1 ) * gShearSubsteps << "," << tDrive[ i ] * tDcrit << "," << tShearForce << ","
                     << tNormalForce << "," << tMeasuredSlip << "," << tShear[ i ] * tLimit << "," << tSlip[ i ] * tDcrit << "\n";
            tLastSlip = tMeasuredSlip;
        }
    }

    void run_case( uint aElements, uint aSubsteps )
    {
        gShearElements      = aElements;
        gShearSubsteps      = aSubsteps;
        gShearFriction      = 0.3;
        gShearYoungsModulus = 1000.0;
        gShearCompression   = 1.0e-4;
        gShearPenalty       = 500.0;
        REQUIRE( par_size() == 1 );
        // Shared EXA main starts console logging. Enable the ASCII trace once,
        // and inspect only this case's portion when Catch runs several cases.
        static bool tLogInitialized = false;
        if ( !tLogInitialized )
        {
            gLogger.initialize( "Coulomb_Shear_solver", 2 );
            tLogInitialized = true;
        }
        std::streampos tLogStart;
        {
            std::ifstream tLog( "Coulomb_Shear_solver.0", std::ios::binary | std::ios::ate );
            REQUIRE( tLog.good() );
            tLogStart = tLog.tellg();
        }

        // A failed workflow must not leave an old Exodus file to satisfy checks.
        std::remove( "Coulomb_Shear.exo" );
        std::remove( "Coulomb_Shear_history.csv" );
        const std::string   tLibrary = MORIS_COULOMB_SHEAR_LIBRARY;
        std::vector< char > tArgument( tLibrary.begin(), tLibrary.end() );
        tArgument.push_back( '\0' );
        char  tName[]      = "";
        char* tArguments[] = { tName, tArgument.data() };
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, tArguments ) == 0 );
        // A workflow can return normally after a Newton iteration limit.
        // Require convergence at every accepted step before examining fields.
        std::ifstream tLog( "Coulomb_Shear_solver.0", std::ios::binary );
        tLog.seekg( tLogStart );
        const std::string tTrace( ( std::istreambuf_iterator< char >( tLog ) ), std::istreambuf_iterator< char >() );
        REQUIRE( tTrace.find( "Hard Stop" ) == std::string::npos );
        uint        tConverged = 0;
        std::size_t tPos       = 0;
        while ( ( tPos = tTrace.find( "Number of Iterations (Convergence)", tPos ) ) != std::string::npos )
        {
            ++tConverged;
            ++tPos;
        }
        REQUIRE( tConverged == tStages * gShearSubsteps );
        check_results();
    }
}    // namespace

TEST_CASE( "Coulomb shear stick-slip cycle", "[CoulombShearCycle]" )
{
    // Compress two elastic layers, then load, slide, hold, unload and reverse.
    // Verify analytic force and slip, friction bound and dissipative direction.
    run_case( 8, 1 );
}

TEST_CASE( "Coulomb shear refined load increments", "[CoulombShearSubdividedLoad]" )
{
    // The same physical path with twice as many accepted increments checks
    // rate-independent Coulomb response and absence of old-stress lag.
    run_case( 8, 2 );
}

TEST_CASE( "Coulomb shear refined mesh", "[CoulombShearRefinedMesh]" )
{
    // Resolve the same two-layer problem on a finer mesh and retain the same
    // analytic physical references and tolerances as the baseline case.
    run_case( 16, 1 );
}
