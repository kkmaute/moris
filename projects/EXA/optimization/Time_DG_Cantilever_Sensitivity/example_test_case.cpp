/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#include <catch.hpp>
#include <cmath>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"
#include <cstdio>
#include <string>
#include <vector>
#include "cl_Communication_Tools.hpp"
#include "HDF5_Tools.hpp"

using namespace moris;

extern "C" {
bool gDGNonlinear          = false;
uint gDGTimeSteps          = 6;
real gDGFinalTime          = 0.3;
real gDGBeamHeight         = 0.43;
real gDGDensity            = 1.0;
real gDGAmplitude          = 0.15;
bool gPrintReferenceValues = false;
}

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

namespace
{
    void run( bool aNonlinear )
    {
        gDGNonlinear                = aNonlinear;
        gDGTimeSteps                = 6;
        gDGFinalTime                = 0.3;
        gDGBeamHeight               = 0.43;
        gDGDensity                  = 1.0;
        gDGAmplitude                = 0.15;
        std::string         library = MORIS_DG_SENSITIVITY_LIBRARY;
        std::vector< char > argument( library.begin(), library.end() );
        argument.push_back( '\0' );
        char  name[]      = "";
        char* arguments[] = { name, argument.data() };
        // Each rank removes its own Exodus output before the workflow so a
        // missing write cannot leave a stale field file for the checker.
        std::string output = "Time_DG_Cantilever_Sensitivity.exo";
        if ( par_size() > 1 )
        {
            std::string rank = std::to_string( par_rank() );
            rank.insert( 0, std::to_string( par_size() ).size() - rank.size(), '0' );
            output += "." + std::to_string( par_size() ) + "." + rank;
        }
        std::remove( output.c_str() );
        if ( par_rank() == 0 )
        {
            // An unsuccessful workflow must not leave an old sweep result
            // available to satisfy the gradient checker.
            std::string sweep = "Time_DG_Cantilever_Sensitivity";
            if ( par_size() > 1 )
            {
                sweep += "_" + std::to_string( par_size() ) + ".0";
            }
            sweep += ".hdf5";
            std::remove( sweep.c_str() );
        }
        barrier();
        REQUIRE( fn_WRK_Workflow_Main_Interface( 2, arguments ) == 0 );
    }

    void check_sensitivities()
    {
        if ( par_rank() != 0 )
        {
            return;
        }
        hid_t            file   = open_hdf5_file( "Time_DG_Cantilever_Sensitivity.hdf5" );
        herr_t           status = 0;
        Matrix< DDRMat > analytical;
        Matrix< DDRMat > volume;
        load_matrix_from_hdf5_file( file, "objective_gradients eval_1-1 analytical", analytical, status );
        REQUIRE( status == 0 );
        load_matrix_from_hdf5_file( file, "constraint_gradients eval_1-1 analytical", volume, status );
        REQUIRE( status == 0 );
        REQUIRE( analytical.numel() == 2 );
        REQUIRE( volume.numel() == 2 );
        // Identify the geometric ADV by dVolume/dh=3. Density has zero
        // volume derivative, independently of GEN's internal ADV ordering.
        const uint height = std::abs( volume( 0 ) ) > std::abs( volume( 1 ) ) ? 0 : 1;
        CHECK( volume( height ) == Approx( 3.0 ).epsilon( 2e-4 ) );
        CHECK( volume( 1 - height ) == Approx( 0.0 ).margin( 1e-9 ) );
        for ( uint p = 0; p < 2; ++p )
        {
            CHECK( std::isfinite( analytical( p ) ) );
            CHECK( std::abs( analytical( p ) ) > 1e-9 );
        }
        for ( uint eps = 1; eps <= 2; ++eps )
        {
            for ( const std::string scheme : { "fd_forward", "fd_backward", "fd_central" } )
            {
                Matrix< DDRMat >  fd;
                const std::string suffix = " eval_1-1 epsilon_" + std::to_string( eps ) + "-2" + " " + scheme;
                load_matrix_from_hdf5_file( file, "objective_gradients" + suffix, fd, status );
                REQUIRE( status == 0 );
                REQUIRE( fd.numel() == 2 );
                for ( uint p = 0; p < 2; ++p )
                {
                    CHECK( std::isfinite( fd( p ) ) );
                    CHECK( analytical( p ) == Approx( fd( p ) ).epsilon( 0.01 ).margin( 1e-7 ) );
                }
            }
        }
        close_hdf5_file( file );
    }
}    // namespace

// Transient linear cantilever under a smooth body-force pulse; compares the
// multislab adjoint against FD for level-set height and material density ADVs.
TEST_CASE( "Time DG linear cantilever sensitivities", "[DGSensitivityLinear]" )
{
    run( false );
    check_sensitivities();
}

// The same geometry, pulse and tolerances with StVK elasticity; checks both
// geometric/material gradients and the velocity-history adjoint recurrence.
TEST_CASE( "Time DG nonlinear cantilever sensitivities", "[DGSensitivityNonlinear]" )
{
    run( true );
    check_sensitivities();
}
