//
// example specific interface to moris
//

#include <catch.hpp>

#include "cl_Logger.hpp" // MRS/IOS/src

//---------------------------------------------------------------

// global variable for interpolation order
uint gInterpolationOrder;

// flag to print reference values
bool gPrintReferenceValues = false;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------

TEST_CASE("Channel_2D_Compressible",
        "[moris],[example],[fluid],[compressible],[Channel_2D_Compressible]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_2D_Compressible.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    REQUIRE( tRet ==  0 );
}
