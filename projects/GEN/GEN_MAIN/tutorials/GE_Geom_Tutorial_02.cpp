/*
 * GO_Geom_Tutorial_02.cpp
 *
 *  Created on: May 21, 2019
 *      Author: sonne
 */

//------------------------------------------------------------------------------
// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
//------------------------------------------------------------------------------
// necessary includes for tutorial
#include "catch.hpp"
// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Element.hpp"
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"

// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

//------------------------------------------------------------------------------

moris::Comm_Manager gMorisComm;

//------------------------------------------------------------------------------
/*!
 * in this tutorial we will:
 *
 * 1) define a matrix of scalar-node values for the discrete representation
 *
 * 2) create a mesh using the HMR database
 *
 * 3) create and initialize a discrete geometry representation class
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
     * <b> Step 1: define a matrix of scalar-node values </b>
     */
    /*!
     * \code{.cpp}
     *     Matrix< DDRMat > tLevelSetVals( 25,1 );
     *     tLevelSetVals(0) =  1.64; tLevelSetVals(1) =   0.64;
     *     tLevelSetVals(2) = -0.36; tLevelSetVals(3) =   0.64;
     *     tLevelSetVals(4) =  0.89; tLevelSetVals(5) =  -0.11;
     *     tLevelSetVals(6) = -0.11; tLevelSetVals(7) =   0.89;
     *     tLevelSetVals(8) =  0.14; tLevelSetVals(9) =   1.64;
     *     tLevelSetVals(10) = 0.64; tLevelSetVals(11) =  0.89;
     *     tLevelSetVals(12) = 0.89; tLevelSetVals(13) = -0.11;
     *     tLevelSetVals(14) = 0.14; tLevelSetVals(15) =  0.64;
     *     tLevelSetVals(16) = 1.64; tLevelSetVals(17) = -0.11;
     *     tLevelSetVals(18) = 0.89; tLevelSetVals(19) =  0.89;
     *     tLevelSetVals(20) = 0.14; tLevelSetVals(21) =  1.64;
     *     tLevelSetVals(22) = 0.89; tLevelSetVals(23) =  0.89;
     *     tLevelSetVals(24) = 0.14;
     * \endcode
     */

    // discrete scalar level set values
    Matrix< DDRMat > tLevelSetVals( 25,1 );
    tLevelSetVals(0) =  1.64; tLevelSetVals(1) =   0.64;
    tLevelSetVals(2) = -0.36; tLevelSetVals(3) =   0.64;
    tLevelSetVals(4) =  0.89; tLevelSetVals(5) =  -0.11;
    tLevelSetVals(6) = -0.11; tLevelSetVals(7) =   0.89;
    tLevelSetVals(8) =  0.14; tLevelSetVals(9) =   1.64;
    tLevelSetVals(10) = 0.64; tLevelSetVals(11) =  0.89;
    tLevelSetVals(12) = 0.89; tLevelSetVals(13) = -0.11;
    tLevelSetVals(14) = 0.14; tLevelSetVals(15) =  0.64;
    tLevelSetVals(16) = 1.64; tLevelSetVals(17) = -0.11;
    tLevelSetVals(18) = 0.89; tLevelSetVals(19) =  0.89;
    tLevelSetVals(20) = 0.14; tLevelSetVals(21) =  1.64;
    tLevelSetVals(22) = 0.89; tLevelSetVals(23) =  0.89;
    tLevelSetVals(24) = 0.14;
//------------------------------------------------------------
    /*!
     * <b> Step 2: create the mesh </b>
     */

    /*!
     * A simple mesh is created from the HMR database. For more information on setting up the mesh, please see the relevant tutorials as the method is
     * only briefly shown here for completeness.
     *
     * \code{.cpp}
     * uint tMeshOrder = 2;
     * \endcode
     *
     * Legrange order = B-spline order = Mesh order = 2
     *
     * \code{.cpp}
     * hmr::ParameterList tParameters = hmr::create_hmr_parameter_list();
     * tParameters.set( "number_of_elements_per_dimension", "2, 2" );
     * tParameters.set( "bspline_orders", "2" );
     * tParameters.set( "lagrange_orders", "2" );
     * tParameters.set( "verbose", 0 );
     * tParameters.set( "domain_dimensions", "2, 2" );
     * tParameters.set( "domain_offset", "-1, -1" );
     * \endcode
     */
    uint tMeshOrder = 2;    // Legrange order = B-spline order = 2

    moris::ParameterList tParameters = hmr::create_hmr_parameter_list();
    tParameters.set( "number_of_elements_per_dimension", "2, 2" );
    tParameters.set( "bspline_orders", "2" );
    tParameters.set( "lagrange_orders", "2" );
    tParameters.set( "verbose", 0 );                    // displays debug info when turned on
    tParameters.set( "domain_dimensions", "2, 2" );     // if no dimensions set, hmr assumes 1x1
    tParameters.set( "domain_offset", "-1, -1" );       // offset so that the center node is at ( 0, 0 )

    /*!
     * Form the parameters set above, the mesh can now be created.
     *
     * \code{.cpp}
     * hmr::HMR tHMR( tParameters );
     * std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tMeshOrder);
     * \endcode
     */
    hmr::HMR tHMR( tParameters );
    // create MTK mesh object and SDF field
    std::shared_ptr< hmr::Mesh > tMesh = tHMR.create_mesh(tMeshOrder);

    /*!
     * Two fields need to be created. One is the field with the scalar values defined above, the other is the field which results from the mapping.
     *
     * \code{.cpp}
     *  std::shared_ptr< hmr::Field > tTargetField = tMesh->create_field( "targetField", tMeshOrder );
     *  std::shared_ptr< hmr::Field > tOutputField = tMesh->create_field( "outputField", tMeshOrder );
     *  \endcode
     */
    std::shared_ptr< hmr::Field > tTargetField = tMesh->create_field( "targetField", tMeshOrder );
    std::shared_ptr< hmr::Field > tOutputField = tMesh->create_field( "outputField", tMeshOrder );

    /*!
     * A mesh pair is then created (assuming the interpolation mesh is the same as the integration mesh) and added to the mesh manager.
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
     * Now the target field is populated with the values defined above and the mesh is finalized.
     *
     * \code{.cpp}
     * tTargetField->put_scalar_values_on_field(tLevelSetVals);
     * tHMR.finalize();
     * \endcode
     */
    tTargetField->put_scalar_values_on_field(tLevelSetVals);
    tHMR.finalize();

    /*!
     * <b> Step 3: create and initialize a discrete representation </b>
     */

    /*!
     * Use the factory to create the geometry representation pointer.
     * \code{.cpp}
     * Ge_Factory tFactory;
     * std::shared_ptr< Geometry > tCircle = tFactory.set_geometry_type(GeomType::DISCRETE);
     * tCircle->set_my_mesh(&tMeshManager);
     * \endcode
     *
     * A discrete representation type will require a target field and an empty output field to be set.
     *
     * \code{.cpp}
     * tCircle->set_my_output_field(tOutputField);
     * tCircle->set_my_target_field(tTargetField);
     * \endcode
     */
    Ge_Factory tFactory;
    std::shared_ptr< Geometry > tCircle = tFactory.set_geometry_type(GeomType::DISCRETE);
    tCircle->set_my_mesh(&tMeshManager);
    tCircle->set_my_output_field(tOutputField);
    tCircle->set_my_target_field(tTargetField);

    /*!
     * <b> Step 4: create the geometry engine and ask for information </b>
     */
    /*!
     * A geometry engine is now created and the discrete representation is set.
     * \code{.cpp}
     * GE_Core tGeomEng;
     * moris_index tMyGeomIndex = tGeomEng.set_geometry( tCircle );
     * \endcode
     *
     * getting the field values are done as follows:
     * \code{.cpp}
     * Matrix< DDRMat > tBSplineCoeff = tGeomEng.get_field_vals(0,15);
     * \endcode
     * This returns the B-Spline coefficient from the first geometry representation set in the GE
     * associated with the 16th B-Spline maxima.
     */
    GE_Core tGeomEng;
    moris_index tMyGeomIndex = tGeomEng.set_geometry( tCircle );

    Matrix< DDRMat > tBSplineCoeff = tGeomEng.get_field_vals(0,15);
    print(tBSplineCoeff, "B-spline coefficient: ");

    //------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;
}






