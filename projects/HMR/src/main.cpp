
// standard
#include <string>

// communication
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MTK_Mesh.hpp"

// core
#include "assert.hpp"
#include "typedefs.hpp"
#include "banner.hpp"

// containers
#include "cl_Cell.hpp"

// linalg

// MTK
#include "cl_MTK_Mesh.hpp"

// HMR
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Fields.hpp"
#include "cl_HMR_State.hpp"
#include "cl_HMR.hpp"

//#include "fn_HMR_EXEC_load_parameters.hpp"

moris::Comm_Manager gMorisComm;

// -----------------------------------------------------------------------------

using namespace moris;
using namespace hmr;

// -----------------------------------------------------------------------------

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
        case( State::INITIALIZE_MESH ) :
        {
            //state_initialize_mesh( tArguments );
            break;
        }
        case( State::REFINE_MESH ) :
        {
            //state_refine_mesh( tArguments );
            break;
        }
        case( State::MAP_FIELDS ) :
        {
            std::cout << "This function is not stable yet" << std::endl;
            //state_map_fields( tArguments );
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
