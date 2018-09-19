
// standard
#include <string>

// communication
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

// core
#include "assert.hpp"
#include "typedefs.hpp"
#include "banner.hpp"

// containers
#include "cl_Cell.hpp"

// linalg
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

// MTK
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Field.hpp"

// STK
#include "cl_STK_Implementation.hpp"

// HMR
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Fields.hpp"
#include "cl_HMR_State.hpp"
#include "cl_HMR.hpp"

moris::Comm_Manager gMorisComm;

// -----------------------------------------------------------------------------

using namespace moris;
using namespace hmr;

// -----------------------------------------------------------------------------

/**
 * This function creates an HMR pointer depending on an XML file and a hdf5 file.
 * Parameters are loaded from the hdf5 file.
 *
 * Parameters except number of elements per direction, domain dimension and origin
 * are overwritten by those set in the XML file.
 *
 */
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

/**
 * this function writes the meshes into specified hdf5 and exo files
 */
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

    // test if exodus outfile is given
    if ( aArguments.get_exodus_output_path().size() > 0 )
    {

        // set exodus timestep
        gStkTimeStep = aArguments.get_timestep();

        // fixme: this functionality is preliminart until STK is fully linked
        // to new MTK
        aHMR->save_to_exodus( aArguments.get_exodus_output_path() );
    }
}

// -----------------------------------------------------------------------------

/**
 * this funciton writes fields into binary files
 */
void
dump_fields(
        Cell< ParameterList >  & aFieldParams,
        Cell< mtk::Field* >    & aOutputFields )
{
    // get number of fields
    uint tNumberOfFields = aFieldParams.size();

    for ( uint f=0; f<tNumberOfFields; ++f )
    {

        // get path to coefficients
        std::string tCoeffPath = aFieldParams( f ).get< std::string >("output_coeffs");

        // get path to node values
        std::string tValuesPath = aFieldParams( f ).get< std::string >("output_values");


        // test if coefficients are to be written
        if ( tCoeffPath.size() > 0 )
        {
            // fixme: this function causes errors in valgrind
            // load values from input
            save_matrix_to_binary_file(
                    aOutputFields( f )->get_coefficients(),
                     tCoeffPath );
        }

        // test if values are to be written
        if ( tValuesPath.size() > 0 )
        {
            // fixme: this function causes errors in valgrind
            // load values from input
            save_matrix_to_binary_file(
                    aOutputFields( f )->get_node_values(),
                    tValuesPath );
        }
    }
}

// -----------------------------------------------------------------------------

/**
 * this function initializes input fields for the refine state routine
 */
void
initialize_input_fields(
              Cell< ParameterList > & aFieldParameters,
              Cell< mtk::Field* >   & aInputFields,
              Mesh                  * aInputMesh,
        const Arguments             & aArguments )
{

    // load field parameters from XML
    load_field_parameters_from_xml(
            aArguments.get_parameter_path(),
            aFieldParameters );

    // get number of fields
    uint tNumberOfFields = aFieldParameters.size();

    // initialize cell of fields
    aInputFields.resize( tNumberOfFields, nullptr );

    // create fields
    for( uint f=0; f<tNumberOfFields; ++f )
    {
        // create field with label
        // fixme: change this if both second and third order are
        //        active at the same time
        aInputFields( f ) = aInputMesh->create_field(
                aFieldParameters( f ).get< std::string >("label") );

        // get path to coefficients
        std::string tCoeffPath = aFieldParameters( f ).get< std::string >("input_coeffs");

        // get path to node values
        std::string tValuesPath = aFieldParameters( f ).get< std::string >("input_values");

        // test if path is set
        if ( tCoeffPath.size() > 0 )
        {
            // load coeffs from file
            load_matrix_from_binary_file(
                    aInputFields( f )->get_coefficients(),
                    tCoeffPath );

            // if values are not loaded from file, calculate them
            if ( tValuesPath.size() == 0 )
            {
                // fixme: HMR needs to calculate T-Matrices for input field first
                MORIS_ERROR( false, "Coefficients can not be loaded yet. HMR will support this feature soon" );

                // calculate node values from coefficients
                aInputFields( f )->evaluate_node_values();
            }
        }
        else if( tValuesPath.size() > 0 )
        {
            // load coeffs from file
            load_matrix_from_binary_file(
                    aInputFields( f )->get_node_values(),
                    tValuesPath );
        }
        else
        {
            MORIS_ERROR( false, "need to provide either path to node values or coefficients to work with input field" );
        }

    }
}

// -----------------------------------------------------------------------------

/**
 * funcition to be called by State::INITIALIZE_MESH
 */
void
state_initialize_mesh( const Arguments & aArguments )
{
    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_hdf5_input_path() );

    // copy input to output
    tHMR->copy_pattern(
            tHMR->get_parameters()->get_input_pattern(),
            tHMR->get_parameters()->get_output_pattern() );

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // delete mesh pointer
    delete tHMR;
}

// -----------------------------------------------------------------------------

/**
 * this funciton is called by State::REFINE_MESH
 */
void
state_refine_mesh( const Arguments & aArguments )
{

    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_hdf5_input_path() );

    tHMR->save_to_exodus( 0, "LastStep.exo" );

    // create pointer to input field
    Mesh * tInputMesh = tHMR->create_input_mesh();

    // cell of input fields
    Cell< mtk::Field* > tInputFields;


    // Cell of parameters for the fields
    Cell< ParameterList > tFieldParameters;

    // initialize fields
    initialize_input_fields(
            tFieldParameters,
            tInputFields,
            tInputMesh,
            aArguments );

    // get number of inout fields
    uint tNumberOfFields = tInputFields.size();

    // counts how many refinemet fields are used
    uint tRefCount = 0;

    // loop over all fields and flag for refinement if flag is set
    for( uint f=0; f<tNumberOfFields; ++f )
    {
        // test if refinement flag is set
        if( tFieldParameters( f ).get< sint >("refine") == 1 )
        {
            // flag volume and surface elements
            tRefCount += tHMR->flag_volume_and_surface_elements( tInputFields( f ) );
        }
    }

    if( tRefCount > 0 ) // can only refine if at least one element was flagged
    {
        // perform refinement routine
        tHMR->perform_refinement();
    }
    else
    {
        // copy input to union
        tHMR->copy_pattern(
                tHMR->get_parameters()->get_input_pattern(),
                tHMR->get_parameters()->get_union_pattern() );

        // copy input to output
        tHMR->copy_pattern(
                tHMR->get_parameters()->get_input_pattern(),
                tHMR->get_parameters()->get_output_pattern() );

        // update meshes
        tHMR->update_meshes();
    }

    // pointer to output mesh
    Mesh * tOutputMesh = tHMR->create_output_mesh();

    // cell of output fields
    Cell< mtk::Field* > tOutputFields( tNumberOfFields, nullptr );

    // project fields to output mesh
    for( uint f=0; f<tNumberOfFields; ++f )
    {
        tOutputFields( f ) = tHMR->map_field_on_mesh( tInputFields( f ), tOutputMesh );
    }

    // write meshes to output
    dump_meshes( aArguments, tHMR );

    // write fields to binary files
    dump_fields( tFieldParameters, tOutputFields );

    // delete input fields
    for( auto tField : tInputFields )
    {
        delete tField;
    }

    // delete interface
    delete tInputMesh;

    // delete fields
    for( auto tField : tOutputFields )
    {
        delete tField;
    }

    // delete interface
    delete tOutputMesh;

    // delete HMR object
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
            state_initialize_mesh( tArguments );
            break;
        }
        case( State::REFINE_MESH ) :
        {
            state_refine_mesh( tArguments );
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
