

#include <string>

#include "cl_Communication_Manager.hpp" // COM/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include "assert.hpp"
#include "cl_Cell.hpp"

#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_State.hpp"
#include "cl_HMR_Parameters.hpp" // HMR/src
#define private public
#define protected public
#include "cl_HMR.hpp" // HMR/src
#undef private
#undef protected

#include "cl_HMR_Fields.hpp"

#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Field.hpp"

#include "fn_load_matrix_from_binary_file.hpp"

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

void
initialize_input_fields(
        const std::string         & aParametersPath,
        const mtk::Mesh           * aMesh,
        Cell< Field_Parameters >  & aFieldParams,
        Cell< mtk::Field* >       & aInputFields )
{
    if( aParametersPath.size() > 0 )
    {

        load_field_parameters_from_xml( aParametersPath, aFieldParams );

        uint tNumberOfFields = aFieldParams.size();

        aInputFields.resize( tNumberOfFields, nullptr );
        for( uint k=0; k<tNumberOfFields; ++k )
        {
            aInputFields( k ) = new mtk::Field(
                    aFieldParams( k ).mLabel,
                    k,
                    aMesh );

            // test if input values are provided
            if ( aFieldParams( k ).mInputValuesPath.size() > 0 )
            {
                // load values from input
                load_matrix_from_binary_file(
                        * aInputFields( k )->get_node_values(),
                        aFieldParams( k ).mInputValuesPath
                        );
            }
        }
    }
}

// -----------------------------------------------------------------------------

void
dump_meshes( const Arguments & aArguments, HMR * aHMR )
{
    // set active pattern to output
    aHMR->set_activation_pattern( aHMR->get_parameters()->get_output_pattern() );

    // test if output path is given
    if ( aArguments.get_hdf5_output_path().size() > 0 )
    {
        aHMR->save_to_hdf5( aArguments.get_hdf5_output_path() );
    }
    //aHMR->mBackgroundMesh->save_to_vtk( "Background.vtk" );

    //aHMR->mLagrangeMeshes( 1 )->save_to_vtk( "Lagrange.vtk" );

    // test if exodus outfile is given
    if ( aArguments.get_exodus_output_path().size() > 0 )
    {
        aHMR->save_to_exodus( aArguments.get_exodus_output_path() );
    }
}

// -----------------------------------------------------------------------------

void
create_tensor_mesh( const Arguments & aArguments )
{
    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_hdf5_input_path() );

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // delete mesh pointer
    delete tHMR;
}

// -----------------------------------------------------------------------------

void
refine_mesh( const Arguments & aArguments )
{
    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_hdf5_input_path() );

    // create pointer to input field
    mtk::Mesh * tInMesh = tHMR->create_mtk_interface( 0 );

    Cell< Field_Parameters >  tFieldParams;
    Cell< mtk::Field* >       tInputFields;

    // initialize fields
    initialize_input_fields(
            aArguments.get_parameter_path(),
            tInMesh,
            tFieldParams,
            tInputFields );

    // get number of fields
    uint tNumberOfFields = tInputFields.size();

    // activate output pattern
    tHMR->set_activation_pattern( tHMR->get_parameters()->get_output_pattern() );

    for( uint k=0; k<tNumberOfFields; ++k )
    {
        if ( tFieldParams( k ).mRefinementFlag )
        {
            Mat<real> tValues = *tInputFields( k )->get_node_values();
            tValues.print("tValues");
            tHMR->flag_against_nodal_field( *tInputFields( k )->get_node_values() );
        }
    }

    tHMR->perform_refinement();

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // delete interface
    delete tInMesh;

    // delete mesh pointer
    delete tHMR;
}

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
            create_tensor_mesh( tArguments );
            break;
        }
        case( State::REFINE_MESH ) :
        {
            refine_mesh( tArguments );
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
