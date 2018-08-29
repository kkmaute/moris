/*
 * cl_XTK_Model_TET10.cpp
 *
 *  Created on: Apr 9, 2018
 *      Author: ktdoble
 */



#include "catch.hpp"


// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/fn_verify_tet_topology.hpp"

// XTKL: Geometry  Include
#include "ios/cl_Logger.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"



#include "xtk/cl_XTK_Model.hpp"

#include "geometry/cl_Sphere.hpp"
#include "geomeng/cl_MGE_Geometry_Engine.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"


namespace xtk
{

TEST_CASE("Generating Tet10s from Tet4s Nonconformal","[TET_10S_NC]")
        {

    // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.50;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 1.0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:2x2x2";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
    tXTKModel.decompose(tDecompositionMethods);

    // verify tet topology before converting to tet10s.
    size_t tCHIndex = 0;
    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tCutMesh = tXTKModel.get_cut_mesh();
    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tChildMesh = tCutMesh.get_child_mesh(0);

    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToNode = tChildMesh.get_element_to_node();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToEdge = tChildMesh.get_element_to_edge();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToFace = tChildMesh.get_element_to_face();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tEdgeToNode    = tChildMesh.get_edge_to_node();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tFaceToNode    = tChildMesh.get_face_to_node();
    bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

    CHECK(tValidTopo);
//    tXTKModel.convert_mesh_tet4_to_tet10();

    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);


    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/xtk_tet10_nonconformal_output.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);

        }


TEST_CASE("Generating Tet10s from Tet4s Conformal","[TET_10S_C]")
        {
    // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);
    tGeometryEngine.mComputeDxDp = false;

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:2x2x1";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    // verify tet topology before converting to tet10s.
    size_t tCHIndex = 0;
    Cut_Mesh<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tCutMesh = tXTKModel.get_cut_mesh();
    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tChildMesh = tCutMesh.get_child_mesh(0);

    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToNode = tChildMesh.get_element_to_node();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToEdge = tChildMesh.get_element_to_edge();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElementToFace = tChildMesh.get_element_to_face();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tEdgeToNode    = tChildMesh.get_edge_to_node();
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tFaceToNode    = tChildMesh.get_face_to_node();
    bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

    CHECK(tValidTopo);

//    tXTKModel.convert_mesh_tet4_to_tet10();


    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);

     std::string tPrefix = std::getenv("XTKOUTPUT");
     std::string tMeshOutputFile = tPrefix + "/xtk_tet10_conformal_output.e";

     tCutMeshData->write_output_mesh(tMeshOutputFile);


        }
}
