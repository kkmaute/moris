/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */


#include "core/cl_XTK_Parameters.hpp"
#include "xtk/cl_XTK_Enums.hpp"

#include "linalg/cl_XTK_Matrix.hpp"

// MPI Header
#include <mpi.h>

// ---------------------------------------------------------------------

namespace xtk
{
Parameters gParameters;
}

int
main( int    argc,
      char * argv[] )
{

    MPI_Init(&argc,&argv);

    xtk::gParameters = xtk::Parameters(Matrix_Backend::EIGEN);


    MPI_Finalize();

    return 0;
}
