/*
 * main.cpp
 *
 *  Created on: Dec 22, 2017
 *      Author: doble
 */

// MORIS header files.

#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Mat.hpp" // LNA/src

moris::Comm_Manager gMorisComm;

using namespace moris;

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // print welcome banner and system information
    moris::print_banner( argc, argv );

//------------------------------------------------------------------------------
//  The main executable is just for developing and testing.
//  Please do not push this file to git.
//------------------------------------------------------------------------------

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
