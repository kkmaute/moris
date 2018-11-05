
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
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

// MTK
#include "cl_MTK_Mesh.hpp"

// HMR
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Field.hpp"
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
            ParameterList tParamList = load_hmr_parameter_list_from_xml( aParametersPath );

            // copy parameters from loaded list ( except offset, number of elements etc )
            aHMR->get_parameters()->copy_selected_parameters( tParamList );
        }

        return aHMR;
    }
    else
    {
        // load parameters from xml path
        ParameterList tParamList = load_hmr_parameter_list_from_xml( aParametersPath );

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
    aHMR->get_database()->set_activation_pattern( aHMR->get_parameters()->get_output_pattern() );

    // test if output path is given
    if ( aArguments.get_database_output_path().size() > 0 )
    {
        aHMR->save_to_hdf5( aArguments.get_database_output_path() );
    }

    // test if exodus outfile is given
    if ( aArguments.get_exodus_output_path().size() > 0 )
    {

        // set exodus timestep

        // to new MTK
        if( aHMR->get_parameters()->get_max_polynomial() < 3 )
        {
            // write debug mesh
            //aHMR->save_mesh_to_vtk("Mesh.vtk");

            // write mesh
            aHMR->save_to_exodus( aArguments.get_exodus_output_path() , aArguments.get_timestep() );


        }
        else
        {

            // activate output pattern
            aHMR->get_database()->set_activation_pattern( aHMR->get_parameters()->get_refined_output_pattern() );

            // write special mesh
            // fixme: the output pattern changes if more that one interpolation is used
            aHMR->save_to_exodus(
                    aHMR->get_parameters()->get_refined_output_pattern(),
                    aArguments.get_exodus_output_path() );
        }
    }

    // write binary coefficients
    if( aArguments.get_binary_path().size() > 0 )
    {
        aHMR->save_coeffs_to_binary_files( aArguments.get_binary_path() );
    }

    // write hdf5 coefficients
    if( aArguments.get_coeffs_path().size() > 0 )
    {
        aHMR->save_coeffs_to_hdf5_file(  aArguments.get_coeffs_path() );
    }
}

// -----------------------------------------------------------------------------

/**
 * this funciton writes fields into binary files
 */
void
dump_fields(
        Cell< ParameterList >            & aFieldParams,
        Cell< std::shared_ptr< Field > > & aFields )
{
    // get number of fields
    uint tNumberOfFields = aFieldParams.size();

    for ( uint f=0; f<tNumberOfFields; ++f )
    {

        // get paths from parameters
        std::string tHDF5FilePath = aFieldParams( f ).get< std::string >("output_hdf5");

        // get path to coefficients
        std::string tCoeffPath    = aFieldParams( f ).get< std::string >("output_coeffs");

        // get path to node values
        std::string tValuesPath   = aFieldParams( f ).get< std::string >("output_values");

        if ( tHDF5FilePath.size() > 0 )
        {
            // load field from HDF5
            aFields( f )->save_field_to_hdf5( tHDF5FilePath );
        }

        if( tCoeffPath.size() > 0 )
        {
            aFields( f )->save_bspline_coeffs_to_binary( tCoeffPath );
        }

        if( tValuesPath.size() > 0 )
        {
            aFields( f )->save_node_values_to_binary( tValuesPath );
        }
    } // end loop over all fields
}

// -----------------------------------------------------------------------------

/**
 * this function initializes input fields for the refine state routine
 */
void
initialize_fields(
        const Arguments                  & aArguments,
        Cell< ParameterList >            & aFieldParameters,
        HMR                              * aHMR,
        Cell< std::shared_ptr< Field > > & aFields )
{

    // load field parameters from XML
    load_field_parameters_from_xml(
            aArguments.get_parameter_path(),
            aFieldParameters );

    // get number of fields
    uint tNumberOfFields = aFieldParameters.size();

    // initialize cell of fields
    aFields.clear();

    // create fields
    for( uint f=0; f<tNumberOfFields; ++f )
    {
        // create field with label
        // fixme: change this if both second and third order are
        //        active at the same time

        // test if input hdf5 is given
        std::string tHDF5FilePath = aFieldParameters( f ).get< std::string >("input_hdf5");

        if ( tHDF5FilePath.size() > 0 )
        {
            sint tLagrange = aFieldParameters( f ).get< sint >( "lagrange_order" );
            sint tBSpline =  aFieldParameters( f ).get< sint >( "bspline_order" );

            // load field from HDF5
            aFields.push_back(aHMR->load_field_from_hdf5_file(
                    tHDF5FilePath,
                    tLagrange,
                    tBSpline ) );
        }
        else
        {

            aFields.push_back(
                    aHMR->create_field(
                            aFieldParameters( f ).get< std::string >("label"),
                            aFieldParameters( f ).get< sint >( "lagrange_order" ),
                            aFieldParameters( f ).get< sint >( "bspline_order" ) ) );

            // get path to coefficients
            std::string tCoeffPath = aFieldParameters( f ).get< std::string >("input_coeffs");

            // get path to node values
            std::string tValuesPath = aFieldParameters( f ).get< std::string >("input_values");

            // test if path is set
            if ( tCoeffPath.size() > 0 )
            {
                // load coeffs from file
                load_matrix_from_binary_file(
                        aFields( f )->get_coefficients(),
                        tCoeffPath );

                // if values are not loaded from file, calculate them
                if ( tValuesPath.size() == 0 )
                {
                    // fixme: HMR needs to calculate T-Matrices for input field first
                    MORIS_ERROR( false, "Coefficients can not be loaded yet. HMR will support this feature soon" );

                    // calculate node values from coefficients
                    aFields( f )->evaluate_node_values();
                }
                else
                {
                    // load coeffs from file
                    load_matrix_from_binary_file(
                            aFields( f )->get_node_values(),
                            tValuesPath );
                }
            }
            else if( tValuesPath.size() > 0 )
            {
                // load coeffs from file
                load_matrix_from_binary_file(
                        aFields( f )->get_node_values(),
                        tValuesPath );
            }
            else
            {
                MORIS_ERROR( false, "You need to provide either path to node values or coefficients to work with input field" );
            }
        }

        // add refinement parameters to field object
        aFields( f )->set_min_volume_level( aFieldParameters( f ).get< sint >("min_volume_refinement_level") );
        aFields( f )->set_max_volume_level( aFieldParameters( f ).get< sint >("max_volume_refinement_level") );
        aFields( f )->set_min_surface_level( aFieldParameters( f ).get< sint >("min_surface_refinement_level") );
        aFields( f )->set_max_surface_level( aFieldParameters( f ).get< sint >("max_surface_refinement_level") );

    } // end loop over all fields
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
            aArguments.get_database_input_path() );

    // get minumum refinement level
    uint tInitialRefinement = tHMR->get_parameters()->get_minimum_initial_refimenent();

    // copy input to output
    if( tInitialRefinement  == 0 )
    {
        tHMR->get_database()->copy_pattern(
                tHMR->get_parameters()->get_input_pattern(),
                tHMR->get_parameters()->get_output_pattern() );
    }
    else
    {

        for( uint k=0; k<tInitialRefinement; ++k )
        {

            // get number of active elements on mesh
            uint tNumberOfElements = tHMR->get_database()->get_number_of_elements_on_proc();

            // flag all elements
            for( uint e=0; e<tNumberOfElements; ++e )
            {
                tHMR->flag_element( e );
            }

            // refine
            tHMR->perform_refinement();
        }
    }
    // special case for third order
    if( tHMR->get_database()->get_parameters()->get_max_polynomial() > 2 )
    {
        tHMR->get_database()->add_extra_refinement_step_for_exodus();
    }

    // finalize database
    tHMR->get_database()->finalize();

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // delete mesh pointer
    delete tHMR;
}

// -----------------------------------------------------------------------------
/**
 * this function is called by State::MAP_FIELDS
 */
void
state_map_fields( const Arguments & aArguments )
{
    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_database_input_path() );

    // cell containing fields
    Cell< std::shared_ptr< Field > > tFields;

    // Cell of parameters for the fields
    Cell< ParameterList > tFieldParameters;

    // load fields from path
    initialize_fields( aArguments, tFieldParameters, tHMR, tFields );

    // check if the user wishes to save the last step
    if( aArguments.get_last_step_path().size() > 0 )
    {
        tHMR->save_last_step_to_exodus(
                aArguments.get_last_step_path(),
                aArguments.get_timestep() );
    }

    // load output database
    tHMR->load_output_pattern_from_path(
            aArguments.get_database_output_path() );


    if( tFields.size() > 0 )
    {
        // map fields
        tHMR->perform_refinement_and_map_fields();

        // save fields into output files
        dump_fields( tFieldParameters, tFields );
    }

    if(  aArguments.get_exodus_output_path().size() > 0 )
    {
        // dump mesh into output
        tHMR->save_to_exodus( aArguments.get_exodus_output_path(), aArguments.get_timestep() );
    }

    // delete mesh pointer
    //delete tHMR;
}
// -----------------------------------------------------------------------------

/**
 * this function is called by State::REFINE_MESH
 */
void
state_refine_mesh( const Arguments & aArguments )
{
    // create mesh pointer
    HMR * tHMR =  initialize_mesh(
            aArguments.get_parameter_path(),
            aArguments.get_database_input_path() );

    // Cell of parameters for the fields
    Cell< ParameterList > tFieldParameters;

    // cell containing fields
    Cell< std::shared_ptr< Field > > tFields;

    // load fields from path
    initialize_fields( aArguments, tFieldParameters, tHMR, tFields );

    // todo: maybe use finalize here?

    // check if the user wishes to save the last step
    if( aArguments.get_last_step_path().size() > 0 )
    {
        tHMR->save_last_step_to_exodus(
                aArguments.get_last_step_path(),
                aArguments.get_timestep() );
    }

    // get number of fields
    uint tNumberOfFields = tFields.size();

    // counts how many elements are refined ( used in the future )
    uint tRefCount = 0;

    // loop over all fields and flag for refinement if flag is set
    for( uint f=0; f<tNumberOfFields; ++f )
    {
        // test if refinement flag is set
        if( tFieldParameters( f ).get< sint >("refine") == 1 )
        {
            uint tCount = tHMR->flag_volume_and_surface_elements( tFields( f ) );

            // flag volume and surface elements
            tRefCount += tCount;


         std::cout << "Refine: proc: " << par_rank() << " field: " << f << " elements: "<<  tCount << std::endl;
        }
    }

    // if no element is flagged for this proc, we flag the parents
    // of all active elements on input pattern
    //if( tRefCount == 0 )
    //{
    //    tHMR->get_database()->get_background_mesh()->flag_active_parents(
    //            tHMR->get_parameters()->get_input_pattern() );
    //}

    tHMR->perform_refinement_and_map_fields();


    if( tNumberOfFields > 0 )
    {
        // save fields into output files
        dump_fields( tFieldParameters, tFields );
    }

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // delete mesh pointer
    //delete tHMR;
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
        case( State::MAP_FIELDS ) :
        {
            std::cout << "This funciton is not stable yet" << std::endl;
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
