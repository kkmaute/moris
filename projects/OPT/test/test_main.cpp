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

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Tools.hpp"    // COM/src
#include "cl_Logger.hpp"                 // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

//------------------------------------------------------------------------------------

void
test_pause()
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

//------------------------------------------------------------------------------------

int
main(
        int   argc,
        char* argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm.initialize( &argc, &argv );

    // for debugging in parallel
    //test_pause();

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    // Run Tests
    int result = Catch::Session().run( argc, argv );

    // finalize moris global communication manager
    gMorisComm.finalize();

    return result;
}
