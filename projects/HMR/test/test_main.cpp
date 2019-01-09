#define CATCH_CONFIG_RUNNER
#include <catch.hpp>


// MORIS header files.
#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src

moris::Comm_Manager gMorisComm;

int
main(
        int    argc,
        char * argv[] )
{
    std::cout << "HMR test made" << std::endl;

    // Initialize Moris global communication manager
    gMorisComm.initialize(&argc, &argv);

    // Run Tests
    int result = Catch::Session().run( argc, argv );

    // finalize moris global communication manager
    gMorisComm.finalize();


    return result;

}


