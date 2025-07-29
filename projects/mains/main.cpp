/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main.cpp
 *
 */

#include <iostream>
#include <ios>
#include <limits>

#include "cl_Communication_Manager.hpp"    // COM/src
#include "cl_Logger.hpp"                   // MRS/IOS/src
#include "banner.hpp"                      // COR/src
#include <Kokkos_Core.hpp>

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

void moris_pause( int& argc, char* argv[] )
{
    // go through user arguments and look for flags
    for ( int k = 0; k < argc; ++k )
    {
        // user requests delayed start
        if ( std::string( argv[ k ] ) == "--pause" || std::string( argv[ k ] ) == "-p" )
        {
            // communicate process ids
            pid_t tPId = getpid();

            Matrix< DDSMat > tPIdVec;
            allgather_scalar( tPId, tPIdVec );

            if ( par_rank() == 0 )
            {
                // print process Ids
                for ( int i = 0; i < par_size(); ++i )
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
    barrier();
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char* argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // delay start of execution for attaching debugger to process
    moris_pause( argc, argv );

    // Kokkos::initialize(argc, argv);

    // set severity level 0 - all outputs
    gLogger.initialize( argc, argv );

    // print banner
    moris::print_banner( argc, argv );

    fprintf( stdout, "par_rank() = %d \n", (int)par_rank() );

    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // Kokkos::finalize_all();

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
