/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * main.cpp
 *
 */

#include <string>
#include <memory>

#include "cl_HMR.hpp"
#include "cl_HMR_Arguments.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Paramfile.hpp"
#include "cl_HMR_State.hpp"
#include "fn_HMR_Exec_dump_fields.hpp"
#include "fn_HMR_Exec_dump_meshes.hpp"
#include "fn_HMR_Exec_initialize_fields.hpp"
#include "fn_HMR_Exec_load_parameters.hpp"
#include "fn_HMR_Exec_load_user_library.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// communication
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

#include "cl_MTK_Mesh.hpp"

// core
#include "assert.hpp"
#include "moris_typedefs.hpp"
#include "banner.hpp"

// containers
#include "cl_Vector.hpp"

// linalg

// MTK
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

// HMR
#include "cl_Logger.hpp" // MRS/IOS/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

using namespace moris;
using namespace hmr;

//// -----------------------------------------------------------------------------
//
//void state_initialize_mesh( const Arguments & aArguments )
//{
//    // create parameter file
//    Paramfile tParams(  aArguments.get_parameter_path(), aArguments.get_state() );
//
//    // create new HMR object from parameter list
//    HMR * tHMR = new HMR( tParams.get_parameter_list() );
//
////    // if there is no initial refinement, copy initial tensor mesh to output
////    if( tHMR->get_parameters()->get_initial_refinement()  == 0 && tHMR->get_parameters()->get_additional_lagrange_refinement() == 0 )
////    {
////        // test if max polynomial is 3
////        if ( tHMR->get_parameters()->get_max_polynomial() > 2 )
////        {
////            // activate extra pattern for exodus
//////            tHMR->get_database()->add_extra_refinement_step_for_exodus();
////        }
////    }
////    else
////    {
////        // otherwise, refine all elements n times
////        tHMR->perform_initial_refinement();
////    }
//
//    // finalize database
//    tHMR->finalize();
//
//    // write mesh
//    dump_meshes( aArguments, tParams, tHMR );
//
//    if (tHMR->get_parameters()->get_refinement_interrelation() )
//    {
//        tHMR->save_mesh_relations_to_hdf5_file( "Mesh_Dependencies.hdf5", 0, 0 );
//    }
//
//    // delete HMR object
//    delete tHMR;
//}
//
//// -----------------------------------------------------------------------------
//
//void state_refine_mesh( const Arguments & aArguments )
//{
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 1: Load Parameter Lists
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // create parameter file
//    Paramfile tParams(  aArguments.get_parameter_path(), aArguments.get_state() );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 2: Initialize HMR Object
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // create new HMR object from parameter list
//    HMR * tHMR = new HMR( tParams.get_input_db_path() );
//
//    // copy parameters from parameter list into HMR object
//    tHMR->get_parameters()->copy_selected_parameters( tParams.get_parameter_list() );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 3: Load user defined functions
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // create library
//    Library tLibrary(  tParams.get_library_path() );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 4: Initialize Fields
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // initialize fields
//    Vector< std::shared_ptr< Field > > tInputFields;
//    initialize_fields( aArguments, tParams, tHMR, tInputFields );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 5: Perform refinement
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // call user defined refinement function
//    //tHMR->based_on_field_put_elements_on_queue( tInputFields, 0, 0 );
//
//    // perform refinement
////    tHMR->perform_refinement_based_on_working_pattern( RefinementMode::LAGRANGE_REFINE );    FIXME
////    tHMR->perform_refinement_based_on_working_pattern( RefinementMode::BSPLINE_REFINE );     FIXME
//
//    // finalize mesh
//    tHMR->finalize();
//
//    //tHMR->save_bsplines_to_vtk("path.vtk");
//
//    if( aArguments.map_while_refine() )
//    {
//        // reserve cell for output fields
//        Vector< std::shared_ptr< Field > > tOutputFields;
//
//        // call mapper
//        perform_mapping( aArguments, tParams, tHMR, tInputFields, tOutputFields );
//
//        // write meshes
//        dump_meshes( aArguments, tParams, tHMR );
//
//        // write fields
//        dump_fields( tParams, tOutputFields );
//    }
//    else
//    {
//        //write meshes
//        dump_meshes( aArguments, tParams, tHMR );
//    }
//
//    if (tHMR->get_parameters()->get_refinement_interrelation() )
//    {
//        tHMR->save_mesh_relations_to_hdf5_file( "Mesh_Dependencies.hdf5", 0, 0 );
//    }
//
//    // delete HMR object
//    delete tHMR;
//}
//// -----------------------------------------------------------------------------
//
//void state_map_fields( const Arguments & aArguments )
//{
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 1: Load Parameter Lists
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // create parameter file
//    Paramfile tParams(  aArguments.get_parameter_path(), aArguments.get_state() );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 2: Initialize HMR Object
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // create new HMR object from parameter list
//    HMR * tHMR = new HMR(
//            tParams.get_input_db_path(),
//            tParams.get_output_db_path() );
//
//    // copy parameters from parameter list into HMR object
//    tHMR->get_parameters()->copy_selected_parameters( tParams.get_parameter_list() );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 3: Initialize Fields
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // initialize fields
//    Vector< std::shared_ptr< Field > > tInputFields;
//    initialize_fields( aArguments, tParams, tHMR, tInputFields );
//
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//    // Step 4: Map Fields and Dump Output
//    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
//
//    // reserve cell for output fields
//    Vector< std::shared_ptr< Field > > tOutputFields;
//
//    // call mapper
//    perform_mapping( aArguments, tParams, tHMR, tInputFields, tOutputFields );
//
//    // write meshes
//    dump_meshes( aArguments, tParams, tHMR );
//
//    // write fields
//    dump_fields( tParams, tOutputFields );
//
//    // delete HMR object
//    delete tHMR;
//}
//
//// -----------------------------------------------------------------------------
//
int
main(
        int    argc,
        char * argv[] )
{
//    // initialize MORIS global communication manager
//    gMorisComm = moris::Comm_Manager( &argc, &argv );
//
//    // Severity level 0 - all outputs
//    gLogger.initialize( 1 );
//
//    // create arguments object
//    Arguments tArguments( argc, argv );
//
//    // select runstate
//    switch ( tArguments.get_state() )
//    {
//        case( State::PRINT_USAGE ) :
//        {
//            // print system usage
//            tArguments.print_usage();
//            break;
//        }
//        case( State::PRINT_VERSION ) :
//        {
//            // print welcome banner and system information
//            moris::print_banner( argc, argv );
//            break;
//        }
//        case( State::PRINT_HELP ) :
//        {
//            // print help line and exit
//            tArguments.print_help();
//            break;
//        }
//        case( State::INITIALIZE_MESH ) :
//        {
//            state_initialize_mesh( tArguments );
//            break;
//        }
//        case( State::REFINE_MESH ) :
//        {
//            state_refine_mesh( tArguments );
//            break;
//        }
//        case( State::MAP_FIELDS ) :
//        {
//            state_map_fields( tArguments );
//            break;
//        }
//        default :
//        {
//            // print system usage
//            tArguments.print_usage();
//            break;
//        }
//    }
//
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
    return 0;
}
//
