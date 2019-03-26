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

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg_typedefs.hpp"

#include "fn_generate_element_to_element.hpp"
#include "fn_generate_shared_face_element_graph.hpp"
#include "fn_assemble_boundary_subphase_constraint.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enums.hpp"
#include "cl_XTK_Cut_Mesh.hpp"

#include "cl_Sphere.hpp"

#include "cl_Matrix.hpp"
#include "cl_MTK_Mesh.hpp"
#include "fn_verify_tet_topology.hpp"

#include "tools/fn_approximate.hpp"


namespace xtk
{

TEST_CASE("Generate shared face element pairs","[SHARED_FACE_ELEM_PAIRS]")
{
    if(par_size() == 1 || par_size() == 2)
    {
    //TODO: FIX THIS TEST in Parallel
    // Geometry Engine Setup ---------------------------------------------------------
    // Using a Levelset Sphere as the Geometry
    real tRadius = 0.25;
    real tXCenter = 1.0;
    real tYCenter = 1.0;
    real tZCenter = 1.0;
    Sphere tLevelsetSphere(tRadius, tXCenter, tYCenter, tZCenter);
    Phase_Table tPhaseTable (1, Phase_Table_Structure::EXP_BASE_2);
    Geometry_Engine tGeometryEngine(tLevelsetSphere,tPhaseTable);

    // Create Mesh --------------------------------------------------------------------
    std::string tMeshFileName = "generated:1x1x2";
    moris::mtk::Mesh* tMeshData = moris::mtk::create_mesh( MeshType::STK, tMeshFileName, nullptr );

    // Setup XTK Model ----------------------------------------------------------------
    size_t tModelDimension = 3;
    Model tXTKModel(tModelDimension,tMeshData,tGeometryEngine);

    //Specify decomposition Method and Cut Mesh ---------------------------------------
    Cell<enum Subdivision_Method> tDecompositionMethods = {Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8};
    tXTKModel.decompose(tDecompositionMethods);

    // Access the Cut Mesh-------------------------------------------------------------
    Cut_Mesh const & tCutMesh = tXTKModel.get_cut_mesh();


    // verify topology of tet 4s
    Child_Mesh tChildMesh = tCutMesh.get_child_mesh(0);

    moris::Matrix< moris::IndexMat > const & tElementToNode = tChildMesh.get_element_to_node();
    moris::Matrix< moris::IndexMat > const & tElementToEdge = tChildMesh.get_element_to_edge();
    moris::Matrix< moris::IndexMat > const & tElementToFace = tChildMesh.get_element_to_face();
    moris::Matrix< moris::IndexMat > const & tEdgeToNode    = tChildMesh.get_edge_to_node();
    moris::Matrix< moris::IndexMat > const & tFaceToNode    = tChildMesh.get_face_to_node();
    bool tValidTopo = verify_tet4_topology(tElementToNode,tElementToEdge,tElementToFace,tEdgeToNode,tFaceToNode);

    CHECK(tValidTopo);


    moris_id tSharedFaceId = 6;
    size_t tParentFaceIndex = tMeshData->get_loc_entity_ind_from_entity_glb_id(tSharedFaceId,EntityRank::FACE);

    size_t tMeshIndex0 = 0;
    size_t tMeshIndex1 = 1;
    size_t tDummy = std::numeric_limits<size_t>::max();

    moris::Matrix< moris::IndexMat > tElementPairs = generate_shared_face_element_pairs(tParentFaceIndex, tMeshIndex0, tMeshIndex1, tCutMesh);

    if(par_size() == 1)

    {
        moris::Matrix< moris::IndexMat > tExpElementPairs(
                {{20, 21, 22, 23},
                 {16, 17, 18, 19}});

        CHECK(equal_to(tElementPairs,tExpElementPairs));
    }

    if(par_size() == 2)
    {
        if(par_rank() == 0)
        {

            moris::Matrix< moris::IndexMat > tExpElementPairs(
                    {{20, 21, 22, 23},
                     {16, 17, 18, 19}});
            CHECK(equal_to(tElementPairs,tExpElementPairs));
        }
        else
        {
            moris::Matrix< moris::IndexMat > tExpElementPairs(
                    {{16, 17, 18, 19},
                     {20, 21, 22, 23}});

            CHECK(equal_to(tElementPairs,tExpElementPairs));
        }
    }
    delete tMeshData;
    }

}
}
#endif /* UNIT_TEST_SRC_XTK_CL_XTK_ELEMENT_TIE_CONSTRAINTS_CPP_ */
