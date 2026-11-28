/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_main.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <Kokkos_Core.hpp>

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src
#include "banner.hpp"                      // COR/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;

//---------------------------------------------------------------

int main( int argc, char *argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // set severity level 0 - all outputs
    gLogger.initialize( 2 );

    // explicitly initialize Kokkos if ArborX is enabled
#ifdef MORIS_HAVE_ARBORX
    Kokkos::initialize( argc, argv );
#endif

    // print banner
    moris::print_banner( argc, argv );

    // Run Tests
    int tRet = Catch::Session().run( argc, argv );

    // finalize Kokkos if ArborX is enabled
#ifdef MORIS_HAVE_ARBORX
    Kokkos::finalize();
#endif

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
