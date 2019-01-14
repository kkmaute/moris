/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */


#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include "xtk/cl_XTK_Enums.hpp"

// MPI Header
#include <mpi.h>

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int
main( int    argc,
      char * argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    int result = 0;
    if(moris::par_size() == 1)
    {
        result = Catch::Session().run( argc, argv );
    }

    // finalize moris global communication manager
    gMorisComm.finalize();

    return result;
}


