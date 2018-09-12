
#include <string>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include "HMR_Arguments.hpp" // HMR/src
#include "cl_HMR_Parameters.hpp" // HMR/src
#include "cl_HMR.hpp" // HMR/src

moris::Comm_Manager gMorisComm;

using namespace moris;
using namespace hmr;
// -----------------------------------------------------------------------------

HMR *
initialize_mesh(
        const std::string & aParametersPath,
        const std::string & aInputPath
         )
{
    // test if an Input Path is given
    if ( aInputPath.size() > 0 )
    {
        // load HMR from file
        HMR * aHMR = new HMR( aInputPath );

        // test if parameters are also given
        if ( aParametersPath.size() > 0 )
        {
            // load parameters from xml path
            ParameterList tParamList = load_parameter_list_from_xml( aParametersPath );

            // copy parameters from loaded list ( except offset, number of elements etc )
            aHMR->get_parameters()->copy_selected_parameters( tParamList );
        }

        return aHMR;
    }
    else
    {

        // load parameters from xml path
        ParameterList tParamList = load_parameter_list_from_xml( aParametersPath );

        // create new HMR object from parameter list
        return new HMR( tParamList );
    }
}

// -----------------------------------------------------------------------------
int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // container for parameter path
    std::string tParametersPath;

    // file with input mesh
    std::string tInputPath;

    // file with output mesh
    std::string tOutputPath;

    // file with exodus output
    std::string tExodusPath;

    // process arguments
    hmr::process_arguments(
            argc,
            argv,
            tParametersPath,
            tInputPath,
            tOutputPath,
            tExodusPath );

    if( tParametersPath.size() > 0 || tInputPath.size() > 0 )
    {
        // print welcome banner and system information
        moris::print_banner( argc, argv );

        HMR * tHMR = initialize_mesh( tParametersPath, tInputPath );


        // test if output path is given
        if ( tOutputPath.size() > 0 )
        {
            tHMR->save_to_hdf5( tOutputPath );
        }

        // test if exodus outfile is given
        if ( tExodusPath.size() > 0 )
        {
            tHMR->save_to_exodus( tExodusPath );
        }

        // delete HMR pointer
        delete tHMR;
    }

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
