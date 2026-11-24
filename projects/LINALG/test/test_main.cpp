/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * test_main.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <chrono>
#include <cstdio>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"    // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int main(
        int   argc,
        char *argv[] )
{
    using Clock       = std::chrono::steady_clock;
    const auto tStart = Clock::now();

    // Keep wall measurements separate from the CPU-time benchmark output.
    // This distinguishes MPI startup/finalization waits from test computation.
    gMorisComm.initialize( &argc, &argv );
    const auto tInitialized = Clock::now();
    const int  tRank        = moris::par_rank();

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // Run Tests
    int        result         = Catch::Session().run( argc, argv );
    const auto tTestsFinished = Clock::now();

    // finalize moris global communication manager
    gMorisComm.finalize();
    const auto tFinished = Clock::now();
    auto       tSeconds  = []( auto aBegin, auto aEnd ) {
        return std::chrono::duration< double >( aEnd - aBegin ).count();
    };
    std::fprintf( stderr,
            "LINALG rank %d wall seconds: MPI_Init=%.6f tests=%.6f MPI_Finalize=%.6f total=%.6f\n",
            tRank,
            tSeconds( tStart, tInitialized ),
            tSeconds( tInitialized, tTestsFinished ),
            tSeconds( tTestsFinished, tFinished ),
            tSeconds( tStart, tFinished ) );

    return result;
}
