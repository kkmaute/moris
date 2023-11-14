/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * HMR_Tutorial_3.cpp
 *
 */

#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
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

//------------------------------------------------------------------------------

/*!
 * \section Tutorial_3: Creating Debug Files
 *
 * Debug files are essential for code learning and development.
 * Apart from the exodus files, which are used for calculation,
 * HMR can also output VTK files for debugging purposes.
 * These files do not contain any fields, but reveal additional
 * data that is usually not visible for the user.
 */
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
//
////------------------------------------------------------------------------------
//
//    /*!
//      * <b> Step 1: create a parameter list </b>
//      */
//
//    /*!
//     * We create a parameterlist and build a mesh with 2x2 elements.
//     * \code{.cpp}
//     * ParameterList tParameters = create_hmr_parameter_list();
//     * tParameters.set( "number_of_elements_per_dimension", "2, 2" );
//     * \endcode
//     */
//    ParameterList tParameters = create_hmr_parameter_list();
//    tParameters.set( "number_of_elements_per_dimension", "4, 4" );
//
//    /*!
//     * Moreover, we want the mesh to be of second order, and again,
//     * we want verbose flags on.
//     * \code{.cpp}
//     * tParameters.set( "bspline_orders", "2" );
//     * tParameters.set( "lagrange_orders", "2" );
//     * tParameters.set( "verbose", 1 );
//     * \endcode
//     */
//    tParameters.set( "bspline_orders", "2" );
//    tParameters.set( "lagrange_orders", "2" );
//
//    /*!
//     * if no dimensions are set, HMR assumes that a 1x1 block is to be built,
//     * and its origin is ( 0, 0 ).
//     */
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 2: create Mesh and perform a manual refinement </b>
//     *
//     * Again, we create a mesh from the parameter list.
//     * \code{.cpp}
//     * HMR tHMR ( tParameters );
//     * \endcode
//     */
//    HMR tHMR ( tParameters );
//
//    /*!
//     * In the following lines, we pick the first element of the mesh,
//     * refine the mesh and repeat the procedure three times.
//     * \code{.cpp}
//     * for( uint k=0; k<3; ++k )
//     * {
//     *   tHMR.flag_element( 0 );
//     *   tHMR.perform_refinement(  moris::hmr::RefinementMode::SIMPLE );
//     * }
//     * \endcode
//     */
//    for( uint k=0; k<3; ++k )
//    {
//        tHMR.flag_element( 0 );
//        tHMR.perform_refinement( moris::hmr::RefinementMode::SIMPLE  );
//    }
//
//
//    /*!
//     * When calling perform_refinement_and_map_fields(), HMR assumes
//     * that one is interested in an output mesh, and also calculates node IDs,
//     * T-Matrices, edges and faces. In contrast, the command perform_refinement()
//     * is meant to be used for intermediate refinement steps.
//     * In order to calculate the data needed for the output, the mesh needs to
//     * \code{.cpp}
//     *   tHMR.finalize();
//     * \endcode
//     */
//    tHMR.finalize();
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 3: Write debug files </b>
//     *
//     */
//
//    /*!
//     * The background mesh contains all elements that exist on the mesh,
//     * including padding elements and, if HMR is run in parallel, elements
//     * on the aura.
//     *
//     * \code{.cpp}
//     * tHMR.save_background_mesh_to_vtk( "BackgroundMesh.vtk" );
//     * \endcode
//     *
//     * The debug file contains the following information
//     *
//     * - ELEMENT_ID : an internal ID where each element is assigned an id
//     *                that is unique to its i-j-k position and its refinement level.
//     *
//     * - ELEMENT_LEVEL : the level on which the element lives
//     *
//     * - ELEMENT_MEMORY_INDEX : the position of the element in the internal memory.
//     *                          HMR also stores deactivated and refined elements.
//     *                          These elements are not visible to MTK.
//     *
//     *
//     * - ELEMENT_OWNER : The proc that owns this element. Note that padding elements
//     *                   are not owned by any proc.
//     */
//    tHMR.save_background_mesh_to_vtk( "BackgroundMesh.vtk" );
//
//    /*!
//     * Unlike the background mesh, which only has elements, the Lagrange mesh
//     * also has nodes. This VTK file does not display padding elements,
//     * but if run in parallel, elements within the aura are shown.
//     *
//     * \code{.cpp}
//     * tHMR.save_mesh_to_vtk( "LagrangeMesh.vtk" );
//     * \endcode
//     *
//     * Node and Element IDs are again unique to level and i-j-k position.
//     * The node and element indices shown are consecutive indices of active
//     * entities over the whole domain. The Entity IDs visible to MTK correspond
//     * to this index plus one, since IDs are 1-based.
//     */
//    tHMR.save_mesh_to_vtk( "LagrangeMesh.vtk" );
//
//
//    /*!
//     * This command dumps the B-Splines into a mesh.
//     *
//     * \code{.cpp}
//     *  tHMR.save_bsplines_to_vtk("BSplines.vtk");
//     * \endcode
//     *
//     * Do the following excercise in ParaView :
//     *
//     * - load the file BSplines.vtk, increase the point size and
//     *   apply a treshold to show only B-Splines that have the state 2.
//     *
//     * - load LagrangeMesh.vtk and display it as surface with edges
//     *
//     * Now, you have visualized the positions of all active B-Splines that are visible to this proc.
//     */
//    tHMR.save_bsplines_to_vtk("BSplines.vtk");
//
//    /*!
//     * This debug file visualizes the facets of this mesh.
//     * In 1D, factets are lines, in 2D, facets are faces.
//     *
//     * \code{.cpp}
//     *  tHMR.save_faces_to_vtk( "Faces.vtk" );
//     * \endcode
//     *
//     */
//    tHMR.save_faces_to_vtk( "Faces.vtk" );
//
//    /*!
//     * The following command dumps the generated edges.
//     * We can't call it in this example, since this mesh is only 2D.
//     *
//     * \code{.cpp}
//     * //tHMR.save_edges_to_vtk("Edges.vtk");
//     * \endcode
//     *
//     */
//
//    //tHMR.save_edges_to_vtk("Edges.vtk");
//
////------------------------------------------------------------------------------
//    // finalize MORIS global communication manager
//    gMorisComm.finalize();
//
//    return 0;
//
}

