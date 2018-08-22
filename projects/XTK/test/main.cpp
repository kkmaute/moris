/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */


#define CATCH_CONFIG_RUNNER
#include "catch.hpp"


#include "core/cl_XTK_Parameters.hpp"
#include "xtk/cl_XTK_Enums.hpp"

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

    int result = Catch::Session().run( argc, argv );

    MPI_Finalize();

    return result;
}


