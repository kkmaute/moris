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
#include "cl_Mesh_Enums.hpp"
#include "fn_verify_tet_topology.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"

// XTKL: Container includes
#include "cl_Cell.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"



#include "cl_XTK_Model.hpp"

#include "cl_Sphere.hpp"
#include "cl_MGE_Geometry_Engine.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"


namespace xtk
{

TEST_CASE("Generating Tet10s from Tet4s Nonconformal","[TET_10S_NC]")
                {
    if(par_size() == 1 || par_size() ==2)
    {
        // Geometry Engine Setup -----------------------
        // Using a Levelset Sphere as the Geometry
        real tRadius = 0.50;
        real tXCenter = 1.0;
        real tYCenter = 1.0;
        real tZCenter = 1.0;
        Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

        Geometry_Engine tGeometryEngine(tLevelsetSphere);

        // Create Mesh ---------------------------------
        std::string tMeshFileName = "generated:2x2x2";
        Cell<std::string> tScalarFields(0);
        mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
        std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

        // Setup XTK Model -----------------------------
        size_t tModelDimension = 3;
        Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

        //Specify your decomposition methods and start cutting
        Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
        tXTKModel.decompose(tDecompositionMethods);

        // verify tet topology before converting to tet10s.
        size_t tCHIndex = 0;
        Cut_Mesh tCutMesh = tXTKModel.get_cut_mesh();
        Child_Mesh tChildMesh = tCutMesh.get_child_mesh(0);

        moris::Matrix< moris::DDSTMat > const & tElementToNode = tChildMesh.get_element_to_node();
        moris::Matrix< moris::DDSTMat > const & tElementToEdge = tChildMesh.get_element_to_edge();
        moris::Matrix< moris::DDSTMat > const & tElementToFace = tChildMesh.get_element_to_face();
        moris::Matrix< moris::DDSTMat > const & tEdgeToNode    = tChildMesh.get_edge_to_node();
        moris::Matrix< moris::DDSTMat > const & tFaceToNode    = tChildMesh.get_face_to_node();
        bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

        CHECK(tValidTopo);
        //    tXTKModel.convert_mesh_tet4_to_tet10();

        std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);


        std::string tPrefix = std::getenv("XTKOUTPUT");
        std::string tMeshOutputFile = tPrefix + "/xtk_tet10_nonconformal_output.e";

        tCutMeshData->write_output_mesh(tMeshOutputFile);
    }
                }


TEST_CASE("Generating Tet10s from Tet4s Conformal","[TET_10S_C]")
{
    if(par_size() == 1 )
    {
    // Geometry Engine Setup -----------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 0;
    Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);

    Geometry_Engine tGeometryEngine(tLevelsetSphere);
    tGeometryEngine.mComputeDxDp = false;

    // Create Mesh ---------------------------------
    std::string tMeshFileName = "generated:2x2x1";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, moris::DDRMat, moris::DDSTMat> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, moris::DDRMat, moris::DDSTMat>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

    // Setup XTK Model -----------------------------
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify your decomposition methods and start cutting
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4};
    tXTKModel.decompose(tDecompositionMethods);

    // verify tet topology before converting to tet10s.
    size_t tCHIndex = 0;
    Cut_Mesh tCutMesh = tXTKModel.get_cut_mesh();
    Child_Mesh tChildMesh = tCutMesh.get_child_mesh(0);

    moris::Matrix< moris::DDSTMat > const & tElementToNode = tChildMesh.get_element_to_node();
    moris::Matrix< moris::DDSTMat > const & tElementToEdge = tChildMesh.get_element_to_edge();
    moris::Matrix< moris::DDSTMat > const & tElementToFace = tChildMesh.get_element_to_face();
    moris::Matrix< moris::DDSTMat > const & tEdgeToNode    = tChildMesh.get_edge_to_node();
    moris::Matrix< moris::DDSTMat > const & tFaceToNode    = tChildMesh.get_face_to_node();
    bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

    CHECK(tValidTopo);

    //    tXTKModel.convert_mesh_tet4_to_tet10();


    std::shared_ptr<mesh::Mesh_Data<xtk::real, xtk::size_t,moris::DDRMat, moris::DDSTMat>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);

    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/xtk_tet10_conformal_output.e";

    tCutMeshData->write_output_mesh(tMeshOutputFile);

    }
}
}
