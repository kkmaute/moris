/*
 * main.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */
//------------------------------------------------------------------------------

// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
#include "cl_Cell.hpp"

//------------------------------------------------------------------------------
// from linalg
#include "cl_Matrix.hpp"
#include "fn_norm.hpp"
#include "fn_load_matrix_from_binary_file.hpp"
#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_print.hpp"

//------------------------------------------------------------------------------
// from MTK
#include "cl_HMR_Field.hpp"

//------------------------------------------------------------------------------

// geometry engine
#include <GEN/src/cl_GEN_Geometry_Engine.hpp>

//------------------------------------------------------------------------------
// HMR
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"

//------------------------------------------------------------------------------
// XTK
#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------



// select namespaces
using namespace moris;
using namespace hmr;
using namespace xtk;

//------------------------------------------------------------------------------
// create communicator
moris::Comm_Manager gMorisComm;


int
main(
        int    argc,
        char * argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm.initialize( &argc, &argv );
//
////------------------------------------------------------------------------------
//
//    /*!
//     * <b> Step 1: Create a Parameter list </b>
//     */
//
//    /*!
//     * We create a parameterlist and build a mesh with 4x4x4 elements.
//     * We want it to be 2x2x2 units long, and the origin at (0,0,0)
//     * \code{.cpp}
//     * ParameterList tParameters = create_hmr_parameter_list();
//     * tParameters.set( "number_of_elements_per_dimension", "4, 4, 4" );
//     * tParameters.set( "domain_dimensions",                "2, 2, 2" );
//     * tParameters.set( "domain_offset",                    "0, 0, 0" );
//     * \endcode
//     */
//    ParameterList tParameters = create_hmr_parameter_list();
//    tParameters.set( "number_of_elements_per_dimension", "2, 2, 2" );
//    tParameters.set( "domain_dimensions",                "2, 2, 2" );
//    tParameters.set( "domain_offset",                    "0.0, 0.0, 0.0" );
//
//
//    /*!
//     * We want to refine the mesh four times. We set the surface refinement
//     * level to 4
//     * \code{.cpp}
//     * tParameters.set( "max_surface_refinement_level", 4 );
//     * \endcode
//     */
//    tParameters.set( "max_surface_refinement_level", 4 );
//
//
//    /*!
//     * Make HMR talkative
//     * \code{.cpp}
//     * tParameters.set( "verbose", 1 );
//     * \endcode
//     */
//    tParameters.set( "verbose", 1 );
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
//    auto tField = tMesh->create_field( "Sphere" );
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
//    for( uint k=0; k<3; ++k )
//    {
//        tField->evaluate_scalar_function( SphereFunction );
//        tHMR.flag_surface_elements( tField );
//        tHMR.perform_refinement();
//        tHMR.update_refinement_pattern();
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
//     * \code{.cpp}
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
//------------------------------------------------------------------------------
//
//    real tRadius  = 95.25;
//    real tXCenter = 1.0;
//    real tYCenter = 1.0;
//    real tZCenter = 0.0;
//    Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
//    Phase_Table tPhaseTable (1,  Phase_Table_Structure::EXP_BASE_2);
//    Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);
//
//    // Create Mesh --------------------------------------------------------------------
//    std::string tMeshFileName = "generated:100x100x30";
//    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, NULL );
//
//    // Setup XTK Model ----------------------------------------------------------------
//    size_t tModelDimension = 3;
//    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);
//    tXTKModel.mVerbose = true;
//
//    //Specify decomposition Method and Cut Mesh ---------------------------------------
//    xtk::Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8,
//                                                                  Subdivision_Method::C_HIERARCHY_TET4};
//    tXTKModel.decompose(tDecompositionMethods);
//


//------------------------------------------------------------------------------
    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return 0;

}
