/*
 * cl_XTK_Cut_Mesh_RegSub.cpp
 *
 *  Created on: Feb 23, 2018
 *      Author: ktdoble
 */
#include "catch.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg/cl_XTK_Matrix.hpp"

#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"
#include "xtk/cl_XTK_Output_Options.hpp"
#include "xtk/fn_mesh_flood_fill.hpp"
#include "xtk/fn_generate_element_to_element.hpp"
#include "xtk/cl_XTK_Enrichment.hpp"
#include "xtk/fn_local_child_mesh_flood_fill.hpp"



// DEBUGGING UTILITY INCLUDES
#include "tools/fn_tet_volume.hpp"
#include "geomeng/fn_Triangle_Geometry.hpp" // For surface normals

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Tools.hpp"
#include "mesh/cl_Mesh_Bucket.hpp"
#include "mesh/cl_Mesh_Side_Set_Input.hpp"
#include "mesh/cl_Mesh_Node_Set_Input.hpp"
#include "mesh/fn_verify_tet_topology.hpp"

namespace xtk
{
TEST_CASE("Direct Testing of the NEW regular subdivision","[NEW_REG_SUB_TEMPLATE]")
{
    // Model dimension
    size_t tModelDim = 3;

    // Set up global coordinates
    moris::Mat_New<real,Default_Matrix_Real> tNodeCoords(15,3);
    tNodeCoords(0,0)  = 0.0; tNodeCoords(0,1)  = 0.0; tNodeCoords(0,2)  = 0.0;
    tNodeCoords(1,0)  = 1.0; tNodeCoords(1,1)  = 0.0; tNodeCoords(1,2)  = 0.0;
    tNodeCoords(2,0)  = 0.0; tNodeCoords(2,1)  = 1.0; tNodeCoords(2,2)  = 0.0;
    tNodeCoords(3,0)  = 1.0; tNodeCoords(3,1)  = 1.0; tNodeCoords(3,2)  = 0.0;
    tNodeCoords(4,0)  = 0.0; tNodeCoords(4,1)  = 0.0; tNodeCoords(4,2)  = 1.0;
    tNodeCoords(5,0)  = 1.0; tNodeCoords(5,1)  = 0.0; tNodeCoords(5,2)  = 1.0;
    tNodeCoords(6,0)  = 0.0; tNodeCoords(6,1)  = 1.0; tNodeCoords(6,2)  = 1.0;
    tNodeCoords(7,0)  = 1.0; tNodeCoords(7,1)  = 1.0; tNodeCoords(7,2)  = 1.0;
    tNodeCoords(8,0)  = 0.5; tNodeCoords(8,1)  = 0.0; tNodeCoords(8,2)  = 0.5;
    tNodeCoords(9,0)  = 1.0; tNodeCoords(9,1)  = 0.5; tNodeCoords(9,2)  = 0.5;
    tNodeCoords(10,0) = 0.5; tNodeCoords(10,1) = 1.0; tNodeCoords(10,2) = 0.5;
    tNodeCoords(11,0) = 0.0; tNodeCoords(11,1) = 0.5; tNodeCoords(11,2) = 0.5;
    tNodeCoords(12,0) = 0.5; tNodeCoords(12,1) = 0.5; tNodeCoords(12,2) = 0.0;
    tNodeCoords(13,0) = 0.5; tNodeCoords(13,1) = 0.5; tNodeCoords(13,2) = 1.0;
    tNodeCoords(14,0) = 0.5; tNodeCoords(14,1) = 0.5; tNodeCoords(14,2) = 0.5;

    // Initialize the Node Indices
    moris::Mat_New<size_t,Default_Matrix_Integer> tNodeIndex({{0, 1, 3, 2, 4, 5, 7, 6, 8, 9, 10, 11, 12, 13, 14}});

    // Intialize Node Ids
    moris::Mat_New<size_t,Default_Matrix_Integer> tNodeId({{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}});

    // Initialize the ancestry
    // Setup the parent tet ancestry (this should be 1 to 1)
    moris::Mat_New<size_t,Default_Matrix_Integer> tNodesAncestry({{0}});
    moris::Mat_New<size_t,Default_Matrix_Integer> tParentEdgeInds({{0,1,2,3,4,5,6,7,8,9,10,11}});
    moris::Mat_New<size_t,Default_Matrix_Integer> tParentEdgeRanks(1,12,1);
    moris::Mat_New<size_t,Default_Matrix_Integer> tParentFaceInds({{0,1,2,3,4,5}});
    moris::Mat_New<size_t,Default_Matrix_Integer> tParentFaceRanks(1,6,2);
    moris::Mat_New<size_t,Default_Matrix_Integer> tElementsAncestry({{0}});

    // Initialize Template
    Mesh_Modification_Template<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tRegSubTemplate(tElementsAncestry(0,0),
                                                                                                       0,
                                                                                                       tNodeIndex,
                                                                                                       tParentEdgeInds,
                                                                                                       tParentEdgeRanks,
                                                                                                       tParentFaceInds,
                                                                                                       tParentFaceRanks,
                                                                                                       TemplateType::REGULAR_SUBDIVISION_HEX8);

    // Initialize child mesh with template
    Child_Mesh_Test<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tRegSubChildMesh(tRegSubTemplate);


    // Check the volume
    moris::Mat_New<size_t,Default_Matrix_Integer> const & tElemToNode = tRegSubChildMesh.get_element_to_node();
    real tVolume = compute_volume_for_multiple_tets(tNodeCoords,tElemToNode);
    CHECK(approximate(tVolume,1.0));


    //
    moris::Mat_New<size_t,Default_Matrix_Integer> tElementPhase(1,24,0);

    size_t tMax = std::numeric_limits<size_t>::max();
    size_t tNumPhases = 2;

    moris::Mat_New<size_t,Default_Matrix_Integer> tActiveElements({{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23}});
    moris::Mat_New<size_t,Default_Matrix_Integer> tIncludedElementMarker(1,24,1);

    // Run flood fill Algorithm
    moris::Mat_New<size_t,Default_Matrix_Integer> tElementSubphase = flood_fill( tRegSubChildMesh.get_element_to_element(),
                                                                      tElementPhase,
                                                                      tActiveElements,
                                                                      tIncludedElementMarker,
                                                                      tNumPhases,
                                                                      tMax,
                                                                      true);


    moris::Mat_New<size_t,Default_Matrix_Integer> tExpElementSubphase(1,24,0);
    CHECK(equal_to(tExpElementSubphase,tElementSubphase));


    bool tValidTopo = verify_tet4_topology(tRegSubChildMesh.get_element_to_node(),
                                           tRegSubChildMesh.get_element_to_edge(),
                                           tRegSubChildMesh.get_element_to_face(),
                                           tRegSubChildMesh.get_edge_to_node(),
                                           tRegSubChildMesh.get_face_to_node());

    CHECK(tValidTopo);

}

}


