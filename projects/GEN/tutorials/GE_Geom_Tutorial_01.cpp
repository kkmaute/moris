/*
 * GE_Geom_Tutorial_01.cpp
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
#include "cl_GE_Factory.hpp"
#include "cl_GE_Node.hpp"
// linalg includes
#include "cl_Matrix.hpp"
#include "fn_all_true.hpp"
#include "fn_equal_to.hpp"
#include "linalg_typedefs.hpp"
#include "op_equal_equal.hpp"
// MTK includes
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Mesh_Factory.hpp"
#include "cl_MTK_Mesh.hpp"
#include "cl_MTK_Mesh_Data_Input.hpp"
#include "cl_MTK_Mesh_Tools.hpp"
#include "cl_MTK_Scalar_Field_Info.hpp"
#include "cl_MTK_Vertex.hpp"

//------------------------------------------------------------------------------

using namespace moris;
using namespace ge;

//------------------------------------------------------------------------------

moris::Comm_Manager gMorisComm;

//------------------------------------------------------------------------------

/*!
 * The following problem is defined:
 *
 * @image html ./figures/geomTutorial_01_setup.jpg "Figure 1: problem setup "
 *
 * In this tutorial we will:
 *
 * 1) create a single-element mesh from mtk using the stk database
 *
 * 2) create an analytic geometry class
 *
 * 3) create the geometry engine and compute information with it
 *
 */

//------------------------------------------------------------------------------

int
main( int    argc,
      char * argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    //------------------------------------------------------------------------------
    /*
     * 1) create a single-element mesh
     * 2) create the Geometry Engine
     * 3) determine LS values and sensitivity at nodes
     * 4) determine intersection locations along edges
     *
     *         [3]
     *   (0,1)      (1,1)
     *      x--------x
     *      |        |
     *      O        |
     * [4]  |        |  [2]
     *      |        |
     *      x-----O--x
     *   (0,0)      (1,0)
     *         [1]
     */

    /*!
     * <b> Step 1: create single-element mesh "by hand" using stk database </b>
     */
    /*!
     * For this tutorial, a single element mesh is defined.
     *
     * For more information on setting up a mesh please see
     * the relevant documentation. Here, the steps are defined only briefly.
     */
    /*!
     * \code{.cpp}
     * uint aNumElemTypes = 1;
     * uint aNumDim = 2;
     *
     * Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};
     * Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};
     * Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
     *                             { 1.0, 0.0 },
     *                             { 1.0, 1.0 },
     *                             { 0.0, 1.0 }};
     *
     * Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }};
     *
     * mtk::MtkMeshData tMeshData( aNumElemTypes );
     * tMeshData.CreateAllEdgesAndFaces = true;
     * tMeshData.SpatialDim = & aNumDim;
     * tMeshData.ElemConn(0) = & aElementConnQuad;
     *
     * tMeshData.NodeCoords = & aCoords;
     * tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
     * tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
     *
     * moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
     * std::string tFieldName = "circle";
     * tNodeField1.set_field_name(tFieldName);
     * tNodeField1.set_field_entity_rank(EntityRank::NODE);
     *
     * moris::mtk::MtkFieldsInfo tFieldsInfo;
     * add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
     * tMeshData.FieldsInfo = &tFieldsInfo;
     * \endcode
     */
    uint aNumElemTypes = 1;     // quad
    uint aNumDim = 2;           // specify number of spatial dimensions

    Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};   // specify element connectivity of quad for mesh

    Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};      // specify the local to global element map for quads

    Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
                                { 1.0, 0.0 },
                                { 1.0, 1.0 },
                                { 0.0, 1.0 }};             // Node coordinate matrix

    Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }}; // specify the local to global map
    //------------------------------------------------------------------------------
    // create MORIS mesh using MTK database
    mtk::MtkMeshData tMeshData( aNumElemTypes );
    tMeshData.CreateAllEdgesAndFaces = true;
    tMeshData.SpatialDim = & aNumDim;
    tMeshData.ElemConn(0) = & aElementConnQuad;
    //------------------------------------------------------------------------------
    tMeshData.NodeCoords = & aCoords;
    tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
    tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
    //------------------------------------------------------------------------------
    // declare scalar node field for the circle LS
    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
    std::string tFieldName = "circle";
    tNodeField1.set_field_name(tFieldName);
    tNodeField1.set_field_entity_rank(EntityRank::NODE);

    // initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    // Place the node field into the field info container
    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);

    // declare some supplementary fields
    tMeshData.FieldsInfo = &tFieldsInfo;

    /*!
     * Now, create the mesh pair to be placed into the mesh manager.
     *
     * For this exercise, the interpolation mesh is the same as the integration mesh.
     *
     * \code{.cpp}
     * mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
     * mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);
     *
     * mtk::Mesh_Manager tMeshManager;
     * uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);
     * \endcode
     */
    // create mesh pair
    mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
    mtk::Integration_Mesh*   tIntegMesh1  = create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);

    /*!
     * <b> Step 2: create the geometry representation </b>
     */
    /*!
     * We first define a moris::Cell of reals which the analytic function uses for required constants (these values are the x,y coordinates of the center
     * and the radius for the circle function being used, see cl_GE_Geometry_Library for other analytic functions).
     *
     * \code{.cpp}
     * moris::Cell< real > tCircleInputs(3);
     * tCircleInputs(0) = 0.0;
     * tCircleInputs(1) = 0.0;
     * tCircleInputs(2) = 0.6;
     * \endcode
     */
    // input parameters for the circle LS
    moris::Cell< real > tCircleInputs(3);
    tCircleInputs(0) = 0.0;   // x center
    tCircleInputs(1) = 0.0;   // y center
    tCircleInputs(2) = 0.6;   // radius

    /*!
     * Create and initialize the geometry representation; depending on the type of geometry representation, the initialization will require different steps.
     * Here, we create an analytic function and pull from the list of known functions in the geometry library.
     *
     * \code{.cpp}
     * Ge_Factory tFactory;
     * std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(GeomType::ANALYTIC);
     * tGeom1->set_analytical_function(AnalyticType::CIRCLE);
     * tGeom1->set_analytical_function_dphi_dx(AnalyticType::CIRCLE);
     * \endcode
     */
    Ge_Factory tFactory;
    std::shared_ptr< Geometry > tGeom1 = tFactory.set_geometry_type(GeomType::ANALYTIC);
    tGeom1->set_analytical_function(AnalyticType::CIRCLE);
    tGeom1->set_analytical_function_dphi_dp(AnalyticType::CIRCLE);
    /*!
     * The mesh needs to be set for all geometry representation types, additionally, the analytic type should have its constants set (otherwise, they are defaulted to zeros).
     * \code{.cpp}
     * tGeom1->set_my_mesh( &tMeshManager );
     * tGeom1->set_my_constants( tCircleInputs );
     * \endcode
     */
    tGeom1->set_my_mesh( &tMeshManager );
    tGeom1->set_my_constants( tCircleInputs );

    /*!
     * <b> Step 3: create the geometry engine and register the created geometry type </b>
     */
    /*!
     * Build the geometry engine and set the geometry we just created.
     * \code{.cpp}
     * GE_Core tGeometryEngine;
     * moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );
     * \endcode
     *
     * When using the .set_geometry() function, a second argument can be passed in which tells the GE the index of the meshes to use in the mesh manager.
     * For this example, there is only one mesh pair, so we use the defaulted value. This implementation would be the same (for this specific case) if we set the geometry using the mesh index from the manager:
     * \code{.cpp}
     * moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1,tMeshIndex );
     * \endcode
     */
    GE_Core tGeometryEngine;
    moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );

    /*!
     * When someone form the outside wants information from the GE, there are two possible ways to approach this:
     *
     * (1) ask the geometry engine directly using the index of the geometry representation object:
     *
     * \code{.cpp}
     *      Matrix< DDRMat > tLSVals                = tGeometryEngine.get_field_vals( tMyGeomIndex );
     *      Cell< Matrix< DDRMat > > tSensitivities = tGeometryEngine.get_sensitivity_vals( tMyGeomIndex );
     *
     * \endcode
     * Where asking this returns the values at all the nodes. If you only want the value at a specific node, pass in the index:
     * \code{.cpp}
     *      Matrix< DDRMat > tLSVals                = tGeometryEngine.get_field_vals( tMyGeomIndex, tVertexIndex );
     *      Cell< Matrix< DDRMat > > tSensitivities = tGeometryEngine.get_sensitivity_vals( tMyGeomIndex, tVertexIndex );
     * \endcode
     *
     * (2) ask the geometry engine for the PDV_Info pointer and use the pointer to access/compute all information directly:
     * \code{.cpp}
     *      PDV_Info* tPDVInfo = tGeometryEngine.get_pdv_info_pointer( tMyGeomIndex );
     *
     *      Matrix< DDRMat > tLSVals                = tPDVInfo->get_field_vals(  );
     *      Cell< Matrix< DDRMat > > tSensitivities = tPDVInfo->get_sensitivity_vals(  );
     * \endcode
     * Where you can ask for information pertaining to a specific node again by passing in the index.
     *
     *
     * For outputting purposes, this field (circle LS) can be added to the STK mesh as:
     * \code{.cpp}
     * tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);
     * \endcode
     */
    PDV_Info* tPDVInfo = tGeometryEngine.get_pdv_info_pointer( tMyGeomIndex );

    Matrix< DDRMat > tLSVals                = tPDVInfo->get_field_vals(  );           // phi
    Cell< Matrix< DDRMat > > tSensitivities = tPDVInfo->get_sensitivity_vals(  );     // dphi/dp

    tInterpMesh1->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tLSVals);      // add the determined values as a field on the mesh (for output purposes)

    /*!
     * <b> Step 4: determine intersection location along edges </b>
     */

    /*!
     * The GE will determine if there is an intersection point with a geometry or any field of interest.
     *
     * Build an intersection object. The object can be of many different types; here we make a line type.
     * \code{.cpp}
     *      Intersection_Object_Line tIntersectionObject;
     * \endcode
     */

    Intersection_Object_Line tIntersectionObject;

    /*!
     *  We now specify the end points of the intersection object directly.
     *  Define global position, time if necessary (here is default constant in time), and "field" values.
     *  The 1st edge has global coordinates ranging from x = 0 to x = 1.
     *  \code{.cpp}
     *  Matrix< DDRMat > tGlobalPos = {{0},{1}};
     *  \endcode
     *  The time for this example is constant so we use the default t = 0 to t = 1.
     *  \code{.cpp}
     *  Matrix< DDRMat > tTHat = {{0},{1}};
     *  \endcode
     *  The field values are taken from the previously determined values.
     *  \code{.cpp}
     *  Matrix< DDRMat > tUHat = {{ tLSVals(0) },{ tLSVals(1)}};
     *  \endcode
     *  Now, set the coordinates and parameter points of the intersection object to be available if needing to work with the FEM module.
     *  \code{.cpp}
     *  tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
     *  \endcode
     *
     */
    Matrix< DDRMat > tGlobalPos = {{0},{1}};
    Matrix< DDRMat > tTHat = {{0},{1}};
    Matrix< DDRMat > tUHat = {{ tLSVals(0) },{ tLSVals(1) }};

    tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );

    /*!
     * The PDV_Info object can now be asked to compute the intersection. If there exists an intersection, an index is returned.
     *
     * \code{.CPP}
     *  moris_index tXInd = tPDVInfo->compute_intersection( &tIntersectionObject );
     * \endcode
     * Again, this can be asked directly through the geometry engine by passing in the index of the specific geometry representation.
     * Now, ask for the intersection point.
     * \code{.cpp}
     * Matrix< F31RMat > tIntersectionAlongX = tPDVInfo->get_intersection_point( &tIntersectionObject, tXInd );
     * \endcode
     *
     * A similar structure is followed to compute the intersection sensitivity with respect to the field or the pdv.
     * \code{.cpp}
     * tPDVInfo->compute_intersection_sensitivity( &tIntersectionObject, tXInd );
     *
     * Matrix< DDRMat > tTemp = tPDVInfo->get_intersection_sensitivity( &tIntersectionObject );
     *
     * Matrix< DDRMat > tTemp2 = tPDVInfo->get_dxgamma_dp( &tIntersectionObject );
     * \endcode
     *
     *
     */

    moris_index tXInd = tPDVInfo->compute_intersection( &tIntersectionObject );

    Matrix< F31RMat > tIntersectionAlongX = tPDVInfo->get_intersection_point( &tIntersectionObject, tXInd );

    tPDVInfo->compute_intersection_sensitivity( &tIntersectionObject, tXInd );

    Matrix< DDRMat > tIntersSensitivity = tPDVInfo->get_intersection_sensitivity( &tIntersectionObject );

    Matrix< DDRMat > tPDVSensitivity = tPDVInfo->get_dxgamma_dp( &tIntersectionObject );

    //------------------------------------------------------------------------------
    /*!
     * Clean up.
     * \code{.cpp}
     * delete tInterpMesh1;
     * delete tIntegMesh1;
     * \endcode
     */
    delete tInterpMesh1;
    delete tIntegMesh1;

    //------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;
}
