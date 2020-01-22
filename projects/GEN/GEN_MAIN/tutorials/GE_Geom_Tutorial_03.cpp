/*
 * GE_Geom_Tutorial_03.cpp
 *
 *  Created on: May 21, 2019
 *      Author: sonne
 */
//------------------------------------------------------------------------------
// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"

// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// LINALG includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

//------------------------------------------------------------------------------

moris::Comm_Manager gMorisComm;

//------------------------------------------------------------------------------
/*!
 * in this tutorial we will:
 *
 * 1) create an MTK mesh using the HMR database
 *
 * 2) define the path to the CAD file
 *
 * 3) create a SDF geometry representation type and initialize it
 *
 * 4) create the geometry engine and get information from it
 *
 */
//------------------------------------------------------------------------------

int
main( int    argc,
      char * argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    //------------------------------------------------------------------------------
    /*!
     * <b> Step 1: create an MTK mesh using the HMR database </b>
     */

    /*!
     * The following parameters are defined to be used in creating the mesh.
     *
     * \code{.cpp}
     * hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
     * tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
     * tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
     * tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );
     * tParameters.set( "lagrange_orders",                  "2" );
     * tParameters.set( "verbose", 1 );
     * \endcode
     */
    hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
    tParameters.set( "number_of_elements_per_dimension", "10, 10, 10" );
    tParameters.set( "domain_dimensions",                "5.6, 2.6, 3.4" );
    tParameters.set( "domain_offset",                    "-4.9, 3.25, -1.7" );
    tParameters.set( "lagrange_orders",                  "2" );
    tParameters.set( "verbose", 1 ); // verbosity flag for debugging and other outputs
    /*!
     * Define: LegrangeOrder = B-SplineOrder = MeshOrder.
     * \code{.cpp}
     *  uint tMeshOrder = 2;
     * \endcode
     */
    uint tMeshOrder = 2;

    /*!
     * Create the mesh anf define the field for the SDF.
     *
     * \code{.cpp}
     * hmr::HMR tHMR( tParameters );
     * std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tMeshOrder);
     * std::shared_ptr< hmr::Field > tField = tMesh->create_field( "SDF", tMeshOrder);
     * tHMR.finalize();
     * \endcode
     */
    hmr::HMR tHMR( tParameters );
    // create MTK mesh object and SDF field
    std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tMeshOrder);
    std::shared_ptr< hmr::Field > tField = tMesh->create_field( "SDF", tMeshOrder);
    tHMR.finalize();

    /*!
     * Create the mesh pair and place in the mesh manager.
     *
     * \code{.cpp}
     * std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tMeshOrder, tHMR.mParameters->get_lagrange_output_pattern());
     * std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh(tMeshOrder, tHMR.mParameters->get_lagrange_output_pattern());
     *
     * mtk::Mesh_Manager tMeshManager;
     * uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpolationMesh.get(),tIntegrationMesh.get());
     * \endcode
     */
    std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh = tHMR.create_interpolation_mesh(tMeshOrder, tHMR.mParameters->get_lagrange_output_pattern());
    std::shared_ptr< moris::hmr::Integration_Mesh_HMR >   tIntegrationMesh   = tHMR.create_integration_mesh(tMeshOrder, tHMR.mParameters->get_lagrange_output_pattern());

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpolationMesh.get(),tIntegrationMesh.get());

    /*!
     * <b> Step 2: define the path to the CAD file </b>
     */

    /*!
     * For this exercise, a file "bracket.obj" is used to create a SDF from.
     *
     * \code{.cpp}
     * std::string tObjectPath = "/projects/HMR/tutorials/bracket.obj";
     * tObjectPath = std::getenv("MORISROOT") + tObjectPath;
     * \endcode
     */
    // get path for STL file to load
    std::string tObjectPath = "/projects/HMR/tutorials/bracket.obj";
    tObjectPath = std::getenv("MORISROOT") + tObjectPath;
//------------------------------------------------------------------------------
    /*!
     * <b> Step 3: create a SDF geometry representation type and initialize it </b>
     */
    /*!
     * Create the geometry representation from the factory.
     *
     * \code{.cpp}
     * Ge_Factory tFactory;
     * std::shared_ptr< Geometry > tSDF = tFactory.set_geometry_type(GeomType::SDF);
     * tSDF->set_my_mesh( &tMeshManager );
     * \endcode
     */
    Ge_Factory tFactory;
    std::shared_ptr< Geometry > tSDF = tFactory.set_geometry_type(GeomType::SDF);
    tSDF->set_my_mesh( &tMeshManager );
    /*!
     * The SDF geometry representation type requires the field created above to be set. Once, set, the SDF has to be initialized.
     *
     * \code{.cpp}
     * tSDF->add_hmr_field( tField );
     * tSDF->initialize_sdf( tObjectPath, tInterpolationMesh );
     * \endcode
     */

    tSDF->add_hmr_field( tField );
    tSDF->initialize_sdf( tObjectPath, tInterpolationMesh );

    /*!
     * <b> Step 4: create the geometry engine and get information from it </b>
     */

    /*!
     * Create the geometry engine and set the SDF representation.
     *
     * \code{.cpp}
     * GE_Core tGeomEng;
     * moris_index tGeomIndex = tGeomEng.set_geometry(tSDF);
     * \endcode
     */
    GE_Core tGeomEng;
    moris_index tGeomIndex = tGeomEng.set_geometry(tSDF);

    /*!
     * We can now ask about the information. Note that asking for the field values at this point will simply return the determined SDF values
     * from the initialization step.
     *
     * Matrix< DDRMat > tTempMat = tGeomEng.get_field_vals(tGeomIndex,tNodeIndex);
     *
     */

    //------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;
}
