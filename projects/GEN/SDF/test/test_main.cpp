#define CATCH_CONFIG_RUNNER

#include <catch.hpp>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int
main(
        int    argc,
        char * argv[] )
{
    // Initialize Moris global communication manager
    gMorisComm = moris::Comm_Manager(&argc, &argv);

    // Severity level 0 - all outputs
    gLogger.initialize( 0 );


    // check if path is set
    std::string tMORISROOT = std::getenv("MORISROOT");
    int result;

    if( tMORISROOT.size() > 0 )
    {
        // Run Tests
        result = Catch::Session().run( argc, argv );
    }
    else
    {
        std::cout << "You need to set the $MORISROOT environment variable" << std::endl
                << "in order to run integration tests." << std::endl
                << "It is most likely /home/$(whoami)/codes/moris." << std::endl;
        result = -1;
    }

    // finalize moris global communication manager
    gMorisComm.finalize();

    return result;

}


