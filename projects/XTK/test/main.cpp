/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main.cpp
 *
 */

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "cl_XTK_Enums.hpp"

// MPI Header
#include <mpi.h>

#include <Kokkos_Core.hpp>

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"                   // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

void pause( int& argc, char* argv[] )
{
    // go through user arguments and look for flags
    for ( int k = 0; k < argc; ++k )
    {
        // user requests delayed start
        if ( std::string( argv[ k ] ) == "--pause" || std::string( argv[ k ] ) == "-p" )
        {
            // communicate process ids
            pid_t tPId = getpid();

            moris::Matrix< moris::DDSMat > tPIdVec;
            moris::comm_gather_and_broadcast( tPId, tPIdVec );

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
    }

    // wait until all processors are ready to continue
    moris::barrier();
}

// ---------------------------------------------------------------------

int main( int argc,
        char* argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm.initialize( &argc, &argv );

    // delay start of execution for attaching debugger to process
    pause( argc, argv );

    // Severity level 0 - all outputs
    gLogger.initialize( 2 );

    Kokkos::initialize( argc, argv );

    int result = 0;

    result = Catch::Session().run( argc, argv );

    // finalize moris global communication manager
    gMorisComm.finalize();

    return result;
}
