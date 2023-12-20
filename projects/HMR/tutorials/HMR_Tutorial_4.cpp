/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tutorial_4.cpp
 *
 */

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "moris_typedefs.hpp"
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
//#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"

//------------------------------------------------------------------------------

// select namespaces
using namespace moris;
using namespace hmr;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;
/*!
 * \section Tutorial_4: Using the MTK Interface
 *
 * In this example, we will refine a mesh according to the surface of
 * a sphere. When we are done, we will create an MTK interface object.
 *
 */
real
SphereFunction( const Matrix< DDRMat > & aPoint )
{
    return norm( aPoint ) - 10;
}

int
main(
        int    argc,
        char * argv[] )
{
//    // initialize MORIS global communication manager
//    gMorisComm = moris::Comm_Manager( &argc, &argv );
//
//    // Severity level 0 - all outputs
//    gLogger.initialize( 0 );
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 1: Create a Parameter list </b>
//     */
//
//    /*!
//     * We create a parameter list and build a mesh with 4x4x4 elements.
//     * We want it to be 2x2x2 units long, and the origin at (0,0,0)
//     * \code{.cpp}
//     * ParameterList tParameters = create_hmr_parameter_list();
//     * tParameters.set( "number_of_elements_per_dimension", "4, 4, 4" );
//     * tParameters.set( "domain_dimensions",                "2, 2, 2" );
//     * tParameters.set( "domain_offset",                    "0, 0, 0" );
//     * \endcode
//     */
//    ParameterList tParameters = create_hmr_parameter_list();
//    tParameters.set( "number_of_elements_per_dimension", "4, 4, 4" );
//    tParameters.set( "domain_dimensions",                "2, 2, 2" );
//    tParameters.set( "domain_offset",                    "0, 0, 0" );
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 2: Create HMR object and refine it to the sphere </b>
//     */
//
//    /*!
//     * Create a field with respect to the Parameterlist.
//     * \code{.cpp}
//     * HMR tHMR( tParameters );
//     * \endcode
//     */
//    HMR tHMR( tParameters );
//
//    /*!
//     * The following two lines create an MTK mesh and a field
//     * that is linked to this mesh.
//     * \code{.cpp}
//     * auto tMesh = tHMR.create_mesh();
//     * auto tField = tMesh->create_field( "Sphere" );
//     * \endcode
//     */
//    auto tMesh = tHMR.create_mesh();
//    auto tField = tMesh->create_field( "Sphere", 1 );
//
//
//    /*!
//     * No we evaluate the field, flag all elements on the surface
//     * and perform the refinement.
//     * The command update_refinement_pattern() is required at this
//     * state. It will be removed in the future.
//     * We repeat this process four times.
//     *
//     * \code{.cpp}
//     * for( uint k=0; k<4; ++k )
//     * {
//     *   // evaluate field
//     *   tField->evaluate_scalar_function( SphereFunction );
//     *   tHMR.flag_surface_elements( tField );
//     *   tHMR.perform_refinement();
//     *   tHMR.update_refinement_pattern();
//     * }
//     * \endcode
//     */
//
//    for( uint k=0; k<4; ++k )
//    {
//        tField->evaluate_scalar_function( SphereFunction );
//        tHMR.flag_surface_elements( tField );
//        tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE );
//        tHMR.update_refinement_pattern( 0 );
//    }
//
//    /*!
//     * We call finalize in order to make the T-Matrices, Surfaces, Edges
//     * and IDs available to MTK.
//     *
//     * \code{.cpp}
//     * tHMR.finalize();
//     * \endcode
//     */
//    tHMR.finalize();
//
//    /**
//     * Evaluate the sphere function one final time. After the mesh is finalized.
//    * \code{.cpp}
//     * tField->evaluate_scalar_function( SphereFunction );
//     * \endcode
//     */
//    tField->evaluate_scalar_function( SphereFunction );
//
//    /*!
//     * Let HMR write the sphere as Exodus file.
//     * \code{.cpp}
//     * tHMR.save_to_exodus( "Sphere.exo" );
//     * \endcode
//     */
//    tHMR.save_to_exodus( "Sphere.exo" );
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 3:Using the MTK API</b>
//     *
//     * Now we can use the MTK API to ask for entity connectivities etc.
//     *
//     * some examples:
//     *
//     * \code{.cpp}
//     * std::cout << "Number of Elements on the Mesh :"    << tMesh->get_num_elems()  << std::endl;
//     * std::cout << "Number of Faces on the Mesh    :"    << tMesh->get_num_faces()  << std::endl;
//     * std::cout << "Number of Edges on the Mesh    :"    << tMesh->get_num_edges()  << std::endl;
//     * std::cout << "Number of Nodes on the Mesh    :"    << tMesh->get_num_nodes()  << std::endl;
//     * std::cout << "Number of DOFs on the Mesh     :"     << tMesh->get_num_coeffs() << std::endl;
//     * \endcode
//     */
//
//    std::cout << "Number of Elements on the Mesh :"    << tMesh->get_num_elems()  << std::endl;
//    std::cout << "Number of Faces on the Mesh    :"    << tMesh->get_num_faces()  << std::endl;
//    std::cout << "Number of Edges on the Mesh    :"    << tMesh->get_num_edges()  << std::endl;
//    std::cout << "Number of Nodes on the Mesh    :"    << tMesh->get_num_nodes()  << std::endl;
//    std::cout << "Number of DOFs on the Mesh     :"    << tMesh->get_num_coeffs( tField->get_bspline_order() ) << std::endl;
//
//
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

