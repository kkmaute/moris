/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tutorial_5.cpp
 *
 */

#include <string>

#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"
#include "paths.hpp"
#include "cl_Logger.hpp"
#include "cl_GlobalClock.hpp" // MRS/IOS/src
#include "cl_Tracer.hpp" // MRS/IOS/src

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_SDF_Generator.hpp"

//------------------------------------------------------------------------------
// HMR

//------------------------------------------------------------------------------
// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
/*!
 * \section Tutorial_5: Create an SDF using the MTK interface and the Geometry engine
 *
 */
real
SphereFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 1.2;
}

int
main(
        int    argc,
        char * argv[] )
{
//    // initialize MORIS global communication manager
//
//    {
//        gMorisComm = moris::Comm_Manager( &argc, &argv );
//    }
//
//    // Severity level 0 - all outputs
//    gLogger.initialize( 0 );
//
////------------------------------------------------------------------------------
//
//    ParameterList tParameters = create_hmr_parameter_list();
//
//
//    // settings for teapot
//    //tParameters.set( "number_of_elements_per_dimension", "20, 12, 10" );
//    //tParameters.set( "domain_dimensions",                "20, 12,  10" );
//    //tParameters.set( "domain_offset",                    "-10, -6, -1" );
//    // std::string tObjectPath = "/projects/HMR/tutorials/utah_teapot.obj";
//
//    tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
//    tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
//    tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );
//    std::string tObjectPath = "/projects/HMR/tutorials/bracket.obj";
//
////------------------------------------------------------------------------------
//    // get path for STL file to load
//    tObjectPath = moris::get_base_moris_dir() + tObjectPath;
//
//    // create SDF generator
//    sdf::SDF_Generator tSdfGen( tObjectPath );
//
//    HMR tHMR( tParameters );
//
//    // create MTK mesh object
//    auto tMesh = tHMR.create_mesh();
//
//
//    for( uint k=0; k<3; ++k )
//    {
//       // matrices with surface element IDs
//       Matrix< IndexMat > tSurfaceElements;
//       tSdfGen.raycast( tMesh, tSurfaceElements );
//       // get number of surface elements
//       uint tNumberOfSurfaceElements = tSurfaceElements.length();
//
//       // loop over all elements
//       for( uint e=0; e<tNumberOfSurfaceElements; ++e )
//       {
//           // manually flag element
//           tHMR.flag_element( tSurfaceElements( e ) );
//       }
//
//       // refine
//       tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//    }
//
//    // calculate T-Matrices etc
//    tHMR.finalize();
//
//    // calculate SDF
//    auto tField = tMesh->create_field( "SDF", 1);
//
////------------------------------------------------------------------------------
//
//    tSdfGen.calculate_sdf( tMesh, tField->get_node_values() );
//
//    tHMR.save_to_exodus( "SDF.exo" );
//
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

