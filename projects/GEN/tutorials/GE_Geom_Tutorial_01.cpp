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
     * We first define a moris::Cell of reals which the analytic functions uses for required constants (these values are the x,y coordinates of the center
     * and the radius for the circle function being used, see Geometry_Library for other analytic functions).
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
    tGeom1->set_analytical_function_dphi_dx(AnalyticType::CIRCLE);
    /*!
     * The mesh needs to be set for all geometry representation types, additionally, the analytic type should have its constants set (otherwise, they are defaulted to zeros).
     * \code{.cpp}
     * tGeom1->set_my_mesh(&tMeshManager);
     * tGeom1->set_my_constants(tCircleInputs);
     * \endcode
     */
    tGeom1->set_my_mesh(&tMeshManager);
    tGeom1->set_my_constants(tCircleInputs);

    /*!
     * <b> Step 3: create the geometry engine and ask about information </b>
     */
    /*!
     * Build the geometry engine and set the geometry we just created.
     * \code{.cpp}
     * GE_Core tGeometryEngine;
     * tGeometryEngine.set_geometry( tGeom1 );
     * \endcode
     *
     * When using the .set_geometry() function, a second argument can be passed in which tells the GE the index of the meshes to use in the mesh manager.
     * For this example, there is only one mesh pair, so we use the defaulted value. This implementation would be the same if we set the geometry using the mesh index from the manager:
     * \code{.cpp}
     * moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1,tMeshIndex );
     * \endcode
     */
    GE_Core tGeometryEngine;
    moris_index tMyGeomIndex = tGeometryEngine.set_geometry( tGeom1 );

    /*!
     * Now ask the geometry engine about information relating to any specific geometry representation which has been set.
     *
     * For example, collect all the values of the circle level-set at the nodes and the sensitivities.
     *
     * \code{.cpp}
     * Cell< Matrix< DDRMat > > tLSVals(4);
     * Cell< Matrix< DDRMat > > tSensitivities(4);
     * for(moris_index n=0; n<4; n++)
     * {
     *      tLSVals(n)        = tGeometryEngine.get_field_vals(tMyGeomIndex,n);
     *      tSensitivities(n) = tGeometryEngine.get_sensitivity_vals(tMyGeomIndex,n);
     * }
     * \endcode
     *
     * The arguments in the .get_field_vals() and .get_sensitivity_vals() are indices corresponding to the geometry representation and the vertex, respectively.
     */
    Cell< Matrix< DDRMat > > tLSVals(4);
    Cell< Matrix< DDRMat > > tSensitivities(4);

    for(moris_index n=0; n<4; n++)
    {
        tLSVals(n)        = tGeometryEngine.get_field_vals(tMyGeomIndex,n);
        tSensitivities(n) = tGeometryEngine.get_sensitivity_vals(tMyGeomIndex,n);
    }

    /*!
     * Additional vertices can be added to the geometry engine to store the information on.
     *
     * \code{.cpp}
     * Node tNewNode(0.5,0.5);
     * tNewNode.set_index( 29 );
     *
     * tGeometryEngine.add_vertex_and_value( tNewNode, tMeshIndex );
     * \endcode
     *
     * The information can then be accessed in the same way as any other vertex information.
     * \code{.cpp}
     * Matrix< DDRMat > tNodeVal = tGeometryEngine.get_field_vals( tMyGeomIndex,tNewNode.get_index() );
     * \endcode
     *
     * Where the first argument is the index of the geometry representation and the second is the vertex index.
     */

    Node tNewNode(0.5,0.5);
    tNewNode.set_index( 29 );

    tGeometryEngine.add_vertex_and_value( tNewNode, tMeshIndex );

    Matrix< DDRMat > tNodeVal = tGeometryEngine.get_field_vals( tMyGeomIndex,tNewNode.get_index() );

    /*!
     * <b> Step 4: determine intersection location along edges </b>
     */

    /*!
     * The GE will determine if there is an intersection point with a geometry or any field of interest.
     *
     * Build an intersection object. The object can be of many different types; here we make a line type.
     * \code{.cpp}
     * Intersection_Object_Line tIntersectionObject;
     * \endcode
     */

    Intersection_Object_Line tIntersectionObject;

    /*!
     *  We now specify the end points of the intersection object directly.
     *  Define global position, time if necessary (here is default constant in time), and "field" values.
     *  \code{.cpp}
     *  Matrix< DDRMat > tGlobalPos = {{0},{1}};
     *  Matrix< DDRMat > tTHat = {{0},{1}};
     *  Matrix< DDRMat > tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(1)(0,0) }};
     *
     *  tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );
     *  \endcode
     */
    Matrix< DDRMat > tGlobalPos = {{0},{1}};
    Matrix< DDRMat > tTHat = {{0},{1}};
    Matrix< DDRMat > tUHat = {{ tLSVals(0)(0,0) },{ tLSVals(1)(0,0) }};

    tIntersectionObject.set_coords_and_param_point( tGeom1, tGlobalPos, tTHat, tUHat );

    /*!
     * Now, ask the geometry engine to compute the intersection of the specified geometry representatin and
     * the intersection object.
     *
     * \code{.cpp}
     * tGeometryEngine.compute_intersection( tMyGeomIndex, &tIntersectionObject );
     * \endcode
     *
     * The intersection object can then provide the intersection point (if it exists).
     *
     * \code{.cpp}
     * Matrix< F31RMat > tIntersectionAlongX = tIntersectionObject.get_intersection_point();
     * \endcode
     */
    tGeometryEngine.compute_intersection( tMyGeomIndex, &tIntersectionObject );

    Matrix< F31RMat > tIntersectionAlongX = tIntersectionObject.get_intersection_point();

    //------------------------------------------------------------------------------
    /*!
     * Clean up after yourself.
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
