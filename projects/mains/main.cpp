#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Logger.hpp"                // MRS/IOS/src
#include "banner.hpp"                   // COR/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int main( int argc, char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // set severity level 0 - all outputs
    gLogger.initialize( 2 );

    // print banner
    moris::print_banner( argc, argv );

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
