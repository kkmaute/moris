#define CATCH_CONFIG_RUNNER
#include <catch.hpp>

// MORIS header files.
#ifdef MORIS_HAVE_PARALLEL
#include <mpi.h>
#endif

#include "linalg.hpp"

// ---------------------------------------------------------------------

// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src

moris::Comm_Manager gMorisComm;

int
main(
        int    argc,
        char * argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // Run Tests
    int result = Catch::Session().run( argc, argv );

    // finalize moris global communication manager
    gMorisComm.finalize();


    return 0;

}


