/*
 * main.cpp
 *
 *  Created on: Sep 30, 2018
 *      Author: messe
 */

// MORIS header files.


#include "cl_Stopwatch.hpp"
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src
#include "cl_Matrix.hpp" // LNA/src

#include "cl_SDF_Triangle_Mesh.hpp"

moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace sdf;

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
    sdf::Triangle_Mesh tMesh( "Part_1.obj" );

    return 0;

}
