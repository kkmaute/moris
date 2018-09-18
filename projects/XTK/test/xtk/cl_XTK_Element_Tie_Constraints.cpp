/*tElm
 * cl_XTK_Element_Tie_Constraints.cpp
 *
 *  Created on: Mar 16, 2018
 *      Author: doble
 */

#ifndef UNIT_TEST_SRC_XTK_CL_XTK_ELEMENT_TIE_CONSTRAINTS_CPP_
#define UNIT_TEST_SRC_XTK_CL_XTK_ELEMENT_TIE_CONSTRAINTS_CPP_

#include "catch.hpp"
#include <limits>

// Linear Algebra Includes

#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/fn_generate_shared_face_element_graph.hpp"
#include "xtk/fn_assemble_boundary_subphase_constraint.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/cl_XTK_Model.hpp"
#include "xtk/cl_XTK_Enums.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"

#include "geometry/cl_Sphere.hpp"

#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/fn_verify_tet_topology.hpp"


#include "tools/fn_approximate.hpp"


namespace xtk
{

TEST_CASE("Generate shared face element pairs","[SHARED_FACE_ELEM_PAIRS]")
{
    // Geometry Engine Setup ---------------------------------------------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 1.0;
    Sphere<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    Phase_Table<size_t, Default_Matrix_Integer> tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh --------------------------------------------------------------------
    std::string tMeshFileName = "generated:1x1x2";
    Cell<std::string> tScalarFields(0);
    mesh::Mesh_Builder_Stk<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tMeshBuilder;
    std::shared_ptr<mesh::Mesh_Data<real, size_t, Default_Matrix_Real, Default_Matrix_Integer>> tMeshData = tMeshBuilder.build_mesh_from_string(tMeshFileName, tScalarFields, true);

    // Setup XTK Model ----------------------------------------------------------------
    size_t tModelDimension = 3;
    Model<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify decomposition Method and Cut Mesh ---------------------------------------
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
    tXTKModel.decompose(tDecompositionMethods);

    // Access the Cut Mesh-------------------------------------------------------------
    Cut_Mesh<real,size_t, Default_Matrix_Real, Default_Matrix_Integer> const & tCutMesh = tXTKModel.get_cut_mesh();


    // verify topology of tet 4s
    Child_Mesh_Test<real, size_t, Default_Matrix_Real, Default_Matrix_Integer> tChildMesh = tCutMesh.get_child_mesh(0);

    moris::Matrix< Default_Matrix_Integer > const & tElementToNode = tChildMesh.get_element_to_node();
    moris::Matrix< Default_Matrix_Integer > const & tElementToEdge = tChildMesh.get_element_to_edge();
    moris::Matrix< Default_Matrix_Integer > const & tElementToFace = tChildMesh.get_element_to_face();
    moris::Matrix< Default_Matrix_Integer > const & tEdgeToNode    = tChildMesh.get_edge_to_node();
    moris::Matrix< Default_Matrix_Integer > const & tFaceToNode    = tChildMesh.get_face_to_node();
    bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

    CHECK(tValidTopo);


    std::shared_ptr<mesh::Mesh_Data<real, size_t,Default_Matrix_Real, Default_Matrix_Integer>> tCutMeshData = tXTKModel.get_output_mesh(tMeshBuilder);

    // Write to exodus file
    std::string tPrefix = std::getenv("XTKOUTPUT");
    std::string tMeshOutputFile = tPrefix + "/tie_test.e";
    tCutMeshData->write_output_mesh(tMeshOutputFile);

    size_t tParentFaceIndex = 5;
    size_t tMeshIndex0 = 0;
    size_t tMeshIndex1 = 1;
    size_t tDummy = std::numeric_limits<size_t>::max();

    moris::Matrix< Default_Matrix_Integer > tElementPhase(1,1,tDummy);
    moris::Matrix< Default_Matrix_Integer > tElementIds(1,1,tDummy);
    moris::Matrix< Default_Matrix_Integer > tElementInds(1,1,tDummy);

    moris::Matrix< Default_Matrix_Integer > tElementPairs = generate_shared_face_element_pairs(tParentFaceIndex, tMeshIndex0, tMeshIndex1, tCutMesh);


    moris::Matrix< Default_Matrix_Integer > tExpElementPairs(
    {{20, 21, 22, 23},
     {16, 17, 18, 19}});

    CHECK(equal_to(tElementPairs,tExpElementPairs));



}



}
#endif /* UNIT_TEST_SRC_XTK_CL_XTK_ELEMENT_TIE_CONSTRAINTS_CPP_ */
