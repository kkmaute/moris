

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

#include "cl_HMR_Field.hpp"
#include "cl_HMR_Fields.hpp"

#include "cl_FEM_IWG_L2.hpp"
#include "cl_MDL_Model.hpp"

#include "cl_MTK_Mesh.hpp"

#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

#include "cl_STK_Implementation.hpp"

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

Mat< uint >
create_block_map(
        const Interface   * aMesh )
{
    Mat< uint >  aMap;
    // clean map
    aMap.set_size( 5, 1, MORIS_UINT_MAX );

    // get number of blocks
    uint tNumberOfBlocks = aMesh->get_number_of_blocks();

    // loop over all blocks
    for( uint k=0; k<tNumberOfBlocks; ++k )
    {
        // add entry
        aMap( aMesh->get_block_by_index( k )->get_interpolation_order() ) = k;
    }

    return aMap;
}

// -----------------------------------------------------------------------------

void
initialize_input_fields(
        const std::string         & aParametersPath,
        Interface                 * aInputMesh,
        Cell< Field_Parameters >  & aFieldParams,
        Cell< Field* >            & aInputFields )
{

    load_field_parameters_from_xml( aParametersPath, aFieldParams );

    // create map
    Mat< uint > tBlockMap = create_block_map( aInputMesh );



    uint tNumberOfFields = aFieldParams.size();

    aInputFields.resize( tNumberOfFields, nullptr );


    for( uint k=0; k<tNumberOfFields; ++k )
    {
        // create new field pointer and link to block with specified order
        aInputFields( k ) = new Field(
                aFieldParams( k ).mLabel,
                k,  aInputMesh->get_hmr_block_by_index(
                        tBlockMap( aFieldParams( k ).mInterpolationOrder )  ) );

        // test if input ceofficients are provided

        // fixme: HMR can't calc values from coeffs yet, since T-Matrix calculation
        //        needs to be moved to another spot
        /*if ( aFieldParams( k ).mInputCoeffsPath.size() > 0 )
        {
            load_matrix_from_binary_file(
                    * aInputFields( k )->get_coefficients(),
                    aFieldParams( k ).mInputCoeffsPath
            );
        } */

        // test if input values are provided
        if ( aFieldParams( k ).mInputValuesPath.size() > 0 )
        {
            // fixme: this function causes errors in valgrind
            // load values from input
            load_matrix_from_binary_file(
                    * aInputFields( k )->get_node_values(),
                    aFieldParams( k ).mInputValuesPath
            );
        }
    }
}

// -----------------------------------------------------------------------------

void
map_output_fields(
        const Cell< Field_Parameters >  & aFieldParams,
        HMR                             * aHMR,
        Interface                       * aInputMesh,
        Cell< Field* >                  & aInputFields,
        Interface                       * aOutputMesh,
        Cell< Field* >                  & aOutputFields )
{

    // fixme: careful: patterns and lagrange meshes are not identical
    //        if multiple interpolation orders are used. This needs to be improved.

    // create union of input and output pattern
    aHMR->unite_patterns(
            aHMR->get_parameters()->get_input_pattern(),
            aHMR->get_parameters()->get_output_pattern(),
            aHMR->get_parameters()->get_union_pattern() );

    // create union Mesh
    Interface * tUnionMesh =  aHMR->create_mtk_interface(
            aHMR->get_parameters()->get_union_pattern() );

    // create maps
    Mat< uint > tOutputMap = create_block_map( aOutputMesh );
    Mat< uint > tUnionMap  = create_block_map( tUnionMesh );

    // get number of fields
    uint tNumberOfFields = aInputFields.size();

    // create output fields
    aOutputFields.resize( tNumberOfFields, nullptr );

    aHMR->set_activation_pattern( aHMR->get_parameters()->get_union_pattern() );

    for( uint k=0; k<tNumberOfFields; ++k )
    {
        // create new field pointer and link to block with specified order
        aOutputFields( k ) = new Field(
                aFieldParams( k ).mLabel,
                k, aOutputMesh->get_hmr_block_by_index(
                        tOutputMap( aFieldParams( k ).mInterpolationOrder ) ) );

        // create union field
        auto tUnionField = new Field(
                aFieldParams( k ).mLabel,
                k, tUnionMesh->get_hmr_block_by_index(
                        tUnionMap( aFieldParams( k ).mInterpolationOrder  ) ) );


        // project data onto output mesh
        aHMR->interpolate_field( aInputFields( k ), tUnionField );

        // create IWG object
        moris::fem::IWG_L2 tIWG;

        // create model
        mdl::Model tModel(
                tUnionMesh,
                tIWG,
                * tUnionField->get_node_values(),
                * aOutputFields( k )->get_coefficients() );

        // evaluate output
        aOutputFields( k )->evaluate_node_values();

        // delete union field
        delete tUnionField;
    }


    delete tUnionMesh;
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
        // set exodus timestep
        gStkTimeStep = aArguments.get_timestep();

        aHMR->save_to_exodus( aArguments.get_exodus_output_path() );
    }
}
// -----------------------------------------------------------------------------

void
dump_fields(
        const Cell< Field_Parameters >  & aFieldParams,
        Cell< Field* >                  & aOutputFields )
{
    // get number of fields
    uint tNumberOfFields = aFieldParams.size();

    for ( uint k=0; k<tNumberOfFields; ++k )
    {
        // test if coefficients are to be written
        if ( aFieldParams( k ).mOutputCoeffsPath.size() > 0 )
        {
            // fixme: this function causes errors in valgrind
            // load values from input
            save_matrix_to_binary_file(
                    *aOutputFields( k )->get_coefficients(),
                    aFieldParams( k ).mOutputCoeffsPath );
        }

        // test if values are to be written
        if ( aFieldParams( k ).mOutputValuesPath.size() > 0 )
        {
            // fixme: this function causes errors in valgrind
            // load values from input
            save_matrix_to_binary_file(
                    *aOutputFields( k )->get_node_values(),
                    aFieldParams( k ).mOutputValuesPath );
        }
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
    Interface * tInputMesh = tHMR->create_mtk_interface(
                      tHMR->get_parameters()->get_input_pattern() );

    Cell< Field_Parameters > tFieldParams;
    Cell< Field* >           tInputFields;
    Cell< Field* >           tOutputFields;

    // initialize fields
    initialize_input_fields(
            aArguments.get_parameter_path(),
            tInputMesh,
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
            tHMR->flag_against_nodal_field( *tInputFields( k )->get_node_values() );
        }
    }

    tHMR->perform_refinement();

    // create pointer to output field
    Interface * tOutputMesh = tHMR->create_mtk_interface(
            tHMR->get_parameters()->get_output_pattern() );

    // map fields
    map_output_fields(
            tFieldParams,
            tHMR,
            tInputMesh,
            tInputFields,
            tOutputMesh,
            tOutputFields );

    // dump mesh into output
    dump_meshes( aArguments, tHMR );

    // dump fields into binary files
    dump_fields( tFieldParams, tOutputFields );

    // delete fields

    for( auto tField : tOutputFields )
    {
        delete tField;
    }

    // delete interface
    delete tOutputMesh;

    // delete input fields
    for( auto tField : tInputFields )
    {
       delete tField;
    }

    // delete interface
    delete tInputMesh;

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
