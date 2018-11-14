// communication
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

#include "banner.hpp"

#include "cl_MTK_Mesh.hpp"

#include "cl_MTK_Mapper_State.hpp"
#include "cl_MTK_Mapper_Arguments.hpp"
//#include "fn_MTK_Mapper_Map_Fields.hpp"

moris::Comm_Manager gMorisComm;

// -----------------------------------------------------------------------------

using namespace moris;
using namespace mapper;

int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // create arguments object
    Arguments tArguments( argc, argv );

    // select runstate
    switch ( tArguments.get_state() )
    {
        case( State::PRINT_USAGE ) :
        {
            // print system usage
            tArguments.print_usage();
            break;
        }
        case( State::PRINT_VERSION ) :
        {
            // print welcome banner and system information
            moris::print_banner( argc, argv );
            break;
        }
        case( State::PRINT_HELP ) :
        {
            // print help line and exit
            tArguments.print_help();
            break;
        }
        case( State::MAP_FIELDS ) :
        {
            //map_fields( tArguments );
            break;
        }
        default :
        {
            // print system usage
            tArguments.print_usage();
            break;
        }
    }


    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;
}

// -----------------------------------------------------------------------------
