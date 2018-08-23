/*
 * main.cpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#include <mpi.h>

#include "cl_Communication_Manager.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include "cl_LINALG_Parameters.hpp"
#include "cl_LINALG_Enums.hpp"

using namespace moris;

Linalg_Parameters gParameters;
Comm_Manager      gMorisComm;

int
main( int    argc,
      char * argv[] )
{

    // initialize MORIS global communication manager
    gMorisComm = Comm_Manager(&argc, &argv);

    // Set up linearl algebra parameters
    gParameters = Linalg_Parameters(moris::Backend_Dense_Matrix::EIGEN_DYNAMIC);


    MPI_Finalize();

    return 0;
}
