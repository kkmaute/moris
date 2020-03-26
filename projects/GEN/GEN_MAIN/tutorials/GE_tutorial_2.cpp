//------------------------------------------------------------------------------
// moris core includes
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"
#include "typedefs.hpp"
//------------------------------------------------------------------------------
// necessary includes for tutorial
#include "catch.hpp"

//------------------------------------------------------------------------------
// GE includes
#include "cl_GE_Core.hpp"
#include "cl_GE_Factory.hpp"

// HMR includes
#include "cl_HMR.hpp"
#include "cl_HMR_Field.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Parameters.hpp"
#include "fn_HMR_Exec_perform_mapping.hpp"

// linalg includes
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
 * 1) create a single-element mesh from the STK database
 *
 * 2) discretize a circle function onto the mesh
 *
 * 3) create the geometry engine
 *
 * 4) ask the geometry engine for information via the output object
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
     * <b> Step 1: create single-element mesh from stk database </b>
     */
    /*!
     * For this tutorial, a single element mesh is defined.
     *
     * For more information on setting up a mesh please see
     * the relevant documentation. Here, the steps are defined only briefly.
     */

    /*!
     * Because the geometry engine requires a T-matrix for transformation from the B-spline basis, we pass in this one which corresponds to second-order B-splines.
     */

    /*!
     *
     * \code{.cpp}
     *
     * Matrix< DDRMat > tTMatrix( 4,4, 0.0 );
     * tTMatrix(0,0) = 0.25;
     * tTMatrix(1,1) = 0.25;
     * tTMatrix(2,2) = 0.25;
     * tTMatrix(3,3) = 0.25;
     *
     * \endcode
     */

    Matrix< DDRMat > tTMatrix( 4,4, 0.0 ); //T-matrix to be used for L2 projection
    tTMatrix(0,0) = 0.25;
    tTMatrix(1,1) = 0.25;
    tTMatrix(2,2) = 0.25;
    tTMatrix(3,3) = 0.25;

    /*!
     * The mesh is created "by hand" using the STK database.
     */

    /*!
     * \code{.cpp}
     * uint aNumElemTypes = 1;
     * uint aNumDim = 2;
     * Matrix< IdMat > aElementConnQuad = {{ 1, 2, 3, 4 }};
     * Matrix< IdMat > aElemLocalToGlobalQuad = {{ 1 }};
     * Matrix< DDRMat > aCoords = {{ 0.0, 0.0 },
     *                             { 1.0, 0.0 },
     *                             { 1.0, 1.0 },
     *                             { 0.0, 1.0 }};
     * Matrix< IdMat > aNodeLocalToGlobal = {{ 1, 2, 3, 4 }};
     *
     * mtk::MtkMeshData tMeshData( aNumElemTypes );
     * tMeshData.CreateAllEdgesAndFaces = true;
     * tMeshData.SpatialDim = & aNumDim;
     * tMeshData.ElemConn(0) = & aElementConnQuad;
     * tMeshData.NodeCoords = & aCoords;
     * tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
     * tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
     *
     * moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
     * std::string tFieldName = "circle";
     * tNodeField1.set_field_name(tFieldName);
     * tNodeField1.set_field_entity_rank(EntityRank::NODE);
     *
     * moris::mtk::Scalar_Field_Info<DDRMat> tNodeField2;
     * std::string tFieldName1 = "projectionVals";
     * tNodeField2.set_field_name(tFieldName1);
     * tNodeField2.set_field_entity_rank(EntityRank::BSPLINE_2);
     *
     * moris::mtk::MtkFieldsInfo tFieldsInfo;
     *
     * add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
     * add_field_for_mesh_input(&tNodeField2,tFieldsInfo);
     *
     * tMeshData.FieldsInfo = &tFieldsInfo;
     *
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
    tMeshData.NodeCoords = & aCoords;
    tMeshData.LocaltoGlobalElemMap(0) = & aElemLocalToGlobalQuad;
    tMeshData.LocaltoGlobalNodeMap = & aNodeLocalToGlobal;
    //------------------------------------------------------------------------------
    // declare scalar node field for the circle LS
    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField1;
    std::string tFieldName = "circle";
    tNodeField1.set_field_name(tFieldName);
    tNodeField1.set_field_entity_rank(EntityRank::NODE);

    moris::mtk::Scalar_Field_Info<DDRMat> tNodeField2;
    std::string tFieldName1 = "projectionVals";
    tNodeField2.set_field_name(tFieldName1);
    tNodeField2.set_field_entity_rank(EntityRank::BSPLINE_2);

    // initialize field information container
    moris::mtk::MtkFieldsInfo tFieldsInfo;
    // Place the node field into the field info container
    add_field_for_mesh_input(&tNodeField1,tFieldsInfo);
    add_field_for_mesh_input(&tNodeField2,tFieldsInfo);

    // declare some supplementary fields
    tMeshData.FieldsInfo = &tFieldsInfo;

    /*!
     * A fully defined mesh consists of both an interpolation mesh and an integration mesh. Here the meshes are the same but they can be different.
     * The meshes are then registered with the Mesh Manager and a mesh index is returned. The Mesh Manager can then be asked for relevant information.
     */

    /*!
     * \code{.cpp}
     *
     * mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
     * mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);
     *
     * mtk::Mesh_Manager tMeshManager;
     * uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);
     *
     * \endcode
     *
     */

    //------------------------------------------------------------------------------
    // create mesh pair
    mtk::Interpolation_Mesh* tInterpMesh1 = create_interpolation_mesh( MeshType::STK, tMeshData );
    mtk::Integration_Mesh*   tIntegMesh1  = mtk::create_integration_mesh_from_interpolation_mesh(MeshType::STK,tInterpMesh1);

    // place the pair in mesh manager
    mtk::Mesh_Manager tMeshManager;
    uint tMeshIndex = tMeshManager.register_mesh_pair(tInterpMesh1,tIntegMesh1);
    //------------------------------------------------------------------------------

    /*!
     * <b> Step 2: discretize a circle function onto the mesh </b>
     */

    /*!
     * For this example we discretize an analytic circle function onto the mesh. The circle function requires a list of constants. A library containing a
     * list of analytic functions can be found in the Geometry_Library class.
     */

    /*!
     *
     * \code{.cpp}
     *
     * moris::Cell< real > tInputs(3);
     * tInputs(0) = 0.0;
     * tInputs(1) = 0.0;
     * tInputs(2) = 0.6;
     *
     * \endcode
     */

    // input parameters for circle
    moris::Cell< real > tInputs(3);
    tInputs(0) = 0.0;   // global x center
    tInputs(1) = 0.0;   // global y center
    tInputs(2) = 0.6;   // radius

    /*!
     * We then project the circle values onto the nodes and add the resulting scalar data to the mesh.
     */

    /*!
     * \code{.cpp}
     * uint tNumNodes = tMeshManager.get_interpolation_mesh(tMeshIndex)->get_num_entities(EntityRank::NODE);
     * Matrix< DDRMat > tNodeVals(1,tNumNodes);
     *
     * for(uint i=0; i<tNumNodes; i++)
     * {
     *      tNodeVals(i) = circle_function( tMeshManager.get_interpolation_mesh(tMeshIndex)->get_node_coordinate(i), tInputs );
     * }
     *
     * tMeshManager.get_interpolation_mesh(tMeshIndex)->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tNodeVals);
     *
     * \endcode
     *
     */

    // compute nodal circle values for mesh
    uint tNumNodes = tMeshManager.get_interpolation_mesh(tMeshIndex)->get_num_entities(EntityRank::NODE);
    Matrix< DDRMat > tNodeVals(1,tNumNodes);

    // collect nodal circle values
    for(uint i=0; i<tNumNodes; i++)
    {
        tNodeVals(i) = circle_function( tMeshManager.get_interpolation_mesh(tMeshIndex)->get_node_coordinate(i), tInputs );
    }
    // add nodal circle values to mesh
    tMeshManager.get_interpolation_mesh(tMeshIndex)->add_mesh_field_real_scalar_data_loc_inds(tFieldName, EntityRank::NODE, tNodeVals);

    /*!
     * <b> Step 3: create the geometry engine </b>
     */

    /*!
     * A Geometry Representatino object pointer is created from the geometry factory.
     * Once this is created, relevant data needs to be set depending on its type (analytic, discrete, SDF, etc.).
     */

    /*!
     *
     * \code{.cpp}
     *
     *  Ge_Factory tFactory;
     *  std::shared_ptr< Geometry > circle = tFactory.set_geometry_type(type::DISCRETE);
     *  circle->set_my_constants( tInputs );
     *
     *  Cell<std::string> tFields(1);
     *  tFields(0) = tFieldName;
     *  circle->set_member_variables(& tMeshManager, tFields);
     *
     * \endcode
     *
     */

    Ge_Factory tFactory;
    std::shared_ptr< Geometry > circle = tFactory.set_geometry_type(GeomType::DISCRETE);
    circle->set_my_constants( tInputs );

    Cell<std::string> tFields(1);   // cell of field names
    tFields(0) = tFieldName;

    circle->set_member_variables(& tMeshManager, tFields);

    /*!
     * Now that the geometry representation is fully defined, we can create the Geometry Engine.
     * When creating a geometry engine, a geometry representation along with its associated Mesh Manager and T-matrix must be passed in.
     *
     */

    /*!
     *
     * \code{.cpp}
     *
     *  GE_Core tGeometryEngine( circle, tMeshManager, tTMatrix );
     *
     *  \endcode
     *
     */

    GE_Core tGeometryEngine( circle, tMeshManager, tTMatrix );

    /*!
     * <b> Step 4: ask the geometry engine for information via the output object </b>
     */

    /*!
     * Once created, the geometry engine initializes the output object which is associated with the input geometry representation. An object which contains all
     * the information can then be created via the index of the geometry representation, as is done here. If this is not preffered, the information can be asked for
     * directly from the geometry engine, again via the geometry representation's index.
     */

    /*!
     *
     * \code{.cpp}
     *
     *  Output_Object tInfo = tGeometryEngine.get_output_object_pointer( 0 );
     *
     *  Matrix< DDRMat > tLSVals = tInfo.get_field_vals();
     *  Matrix< DDRMat > tNodalADVs = tInfo.get_nodal_advs();
     *
     *  \endcode
     *
     */

    Output_Object tInfo = tGeometryEngine.get_output_object_pointer( 0 ); // get the object which has all the information

    Matrix< DDRMat > tLSVals = tInfo.get_field_vals();
    Matrix< DDRMat > tNodalADVs = tInfo.get_nodal_advs();


    //------------------------------------------------------------------------------
    gMorisComm.finalize();

    return 0;
}
