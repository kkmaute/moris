/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * query_main.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Logger.hpp"

// MORIS header files.
#include "ios.hpp"

// query function definitions
#include <string>

// Tracer Enums to understand keywords
#include "cl_Tracer_Enums.hpp"

#include "cl_Query.hpp"

// initialize global variables
moris::Logger gLogger;
moris::Comm_Manager gMorisComm;

int
main(
        int    argc,
        char * argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm.initialize(&argc, &argv);

// ----------------------------------------------------- //

        // create query
        moris::ios::Query tQuery;

        // initialize with userinput
        tQuery.run(argc,argv);

// ----------------------------------------------------- //
//    // Running Tests with catch
//    int result = Catch::Session().run( argc, argv );
// ----------------------------------------------------- //

    // finalize moris global communication manager
    gMorisComm.finalize();

    return EXIT_SUCCESS;

}

