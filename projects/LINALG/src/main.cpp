/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main.cpp
 *
 */

#include <mpi.h>

#include "cl_Communication_Manager.hpp" // COM/src
#include "moris_typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include <iostream>
#include <ctime>

#include "linalg_typedefs.hpp"
// Global variables
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;

int
main( int    argc,
      char * argv[] )
{

    // initialize MORIS global communication manager
    gMorisComm = Comm_Manager(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );

    gMorisComm.finalize();

    return 0;
}

