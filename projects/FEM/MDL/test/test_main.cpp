/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 */
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// MORIS header files.
#include "cl_Communication_Manager.hpp"
#include "cl_Logger.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int
main(
        int   argc,
        char* argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm.initialize( &argc, &argv );

    // Severity level 0 - all outputs
    //    gLogger.initialize(0);
    gLogger.initialize( "gLogger.log", 2, 1 );

    // Run Tests
    int result = Catch::Session().run( argc, argv );

    // finalize moris global communication manager
    gMorisComm.finalize();

    return result;
}
