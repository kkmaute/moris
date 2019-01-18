/*
 * cl_XTK_Cut_Mesh_Reg_to_Nh.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "xtk/cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"


// DEBUGGING UTILITY INCLUDES
#include "geomeng/fn_Triangle_Geometry.hpp"


namespace xtk
{


TEST_CASE("Child Mesh","[Child_Mesh]")
{
//
//    size_t tMax = std::numeric_limits<size_t>::max();
//    size_t tParentElementIndex = 1;
//
//    // Processor local node indices
//    moris::Matrix< moris::IndexMat > tNodeInds({{1,2,3,4,5,6}});
//
//    // Locally indexed
//    moris::Matrix< moris::IndexMat > tElementToNode({{1,2,3,4},
//                                                     {2,5,3,4},
//                                                     {1,2,3,6},
//                                                     {1,4,2,5}});
//
//    // Processor local index
//   moris::Matrix< moris::IndexMat > tElementEdgeParentInds({{0, 1, 2, 3, 4, 5},
//                                                              {6, 7, 1, 4, 8, 5},
//                                                              {0, 1, 2, 9,10,11},
//                                                              {6, 7, 1,10,13,14}});
//
//   moris::Matrix< moris::DDSTMat > tElementEdgeParentRanks(4,6,1); // All rank 1
//
//   // Processor local index
//   moris::Matrix< moris::IndexMat > tElementFaceParentInds({{0, 1, 2, 3},
//                                                              {4, 5, 6, 1},
//                                                              {0, 7, 8, 9},
//                                                              {4,10,11, 7}});
//
//   moris::Matrix< moris::DDSTMat > tElementFaceParentRanks(4,4,2); // All rank 2
//   moris::Matrix< moris::DDSTMat > tElementInferfaceSides(4,1,tMax);
//
//
//
//   Child_Mesh_Test tChildMesh(tParentElementIndex,
//                                                                                      tNodeInds,
//                                                                                      tElementToNode,
//                                                                                      tElementEdgeParentInds,
//                                                                                      tElementEdgeParentRanks,
//                                                                                      tElementFaceParentInds,
//                                                                                      tElementFaceParentRanks,
//                                                                                      tElementInferfaceSides);
//
//
//   // BEGIN TEMPLATE ---------------------------------------------------------------------------------------
//   // Insert a template
//   Mesh_Modification_Template<real,size_t,moris::DDRMat,moris::DDSTMat> tMeshTemplate;
//   tMeshTemplate.mNumNewElem = 2;
//   tMeshTemplate.mNumElemToReplace = 1;
//   tMeshTemplate.mElemIndToReplace = 2;
//   tMeshTemplate.mParentEdgeInds           = moris::Matrix< moris::IndexMat >({{2,4,5,7,8,9},
//                                                                               {0,2,1,4,3,9}});
//   tMeshTemplate.mParentEdgeRanks          = moris::Matrix< moris::DDSTMat >(2,6,1);
//   tMeshTemplate.mParentFaceInds           = moris::Matrix< moris::IndexMat >({{2,4,5,7},
//                                                                                 {0,2,1,4}});
//   tMeshTemplate.mParentFaceRanks          = moris::Matrix< moris::DDSTMat >(2,4,2);
//   tMeshTemplate.mNodeInds                 = moris::Matrix< moris::IndexMat >({{10,11,12,13,14,15,16}});
//
//   tMeshTemplate.mNewElementToNode         = moris::Matrix< moris::IndexMat >({{0,1,2,3},
//                                                                                 {4,5,6,1}});
//   tMeshTemplate.mNewParentEdgeRanks       = moris::Matrix< moris::DDSTMat >(2,6,1);
//
//   tMeshTemplate.mNewParentEdgeOrdinals    = moris::Matrix< moris::IndexMat >({{0,1,2,3,4,5},
//                                                                                 {0,1,2,3,4,5}});
//   tMeshTemplate.mNewParentFaceRanks       = moris::Matrix< moris::DDSTMat >(2,4,2);
//
//   tMeshTemplate.mNewParentFaceOrdinals    = moris::Matrix< moris::IndexMat >({{0,1,2,3},
//                                                                                 {3,2,1,0}});
//   tMeshTemplate.mNewElementInterfaceSides = moris::Matrix< moris::DDSTMat >({{1},
//                                                                             {tMax}});
//      // END TEMPLATE -----------------------------------------------------------------------------------------
//
//   tChildMesh.add_node_indices(tMeshTemplate.mNodeInds);
//
//   // Insert template into mesh
//   tChildMesh.allocate_more_elements(tMeshTemplate.mNumNewElem - tMeshTemplate.mNumElemToReplace);
//   tChildMesh.insert_child_mesh_template(tMeshTemplate);
//   tChildMesh.generate_connectivities(true,true,true);
}


/*
 * Test the transition between a regular subdivision hex to conformal node hierarchy mesh independently of the geometry engine
 */
TEST_CASE("Regular Subdivision to NH Transition","[RS_NH]")
{
    moris::Matrix< moris::IndexMat > tRegSubNodes({{0 ,4 ,6 ,2 ,1 ,5 ,7 ,3 ,8 ,9 ,10,11,12,13,14}});
    moris::Matrix< moris::IdMat > tRegSubNodeIds({{1,2,4,3,5,6,8,7,9,10,11,12,13,14,15}});
    moris::Matrix< moris::IdMat > tFullNodeIds({{1,2 ,4 ,3 ,5 ,6 ,8 ,7 ,9 ,10,11,12,13,14,15,19,16,22,17,20,18,21}});

    /*
     * Setup node coordinate list corresponding to the edge indices above
     */
    moris::Matrix< moris::DDRMat > tNodeCoords({{0, 0, 0},
                                                {0, 0, 1},
                                                {0, 1, 0},
                                                {0, 1, 1},
                                                {1, 0, 0},
                                                {1, 0, 1},
                                                {1, 1, 0},
                                                {1, 1, 1},
                                                {0.5, 0, 0.5},
                                                {1, 0.5, 0.5},
                                                {0.5, 1, 0.5},
                                                {0, 0.5, 0.5},
                                                {0.5, 0.5, 0},
                                                {0.5, 0.5, 1},
                                                {0.5, 0.5, 0.5},
                                                {1, 1, 0.6975},
                                                {1, 0.6975, 1},
                                                {0.6975, 1, 1},
                                                {1, 0.6975, 0.6975},
                                                {0.6975, 1, 0.6975},
                                                {0.6975, 0.6975, 1},
                                                {0.7983333333333332, 0.7983333333333332, 0.7983333333333332}});


    /*
     * Inheritance
     */
     moris::Matrix< moris::IndexMat > tParentNodeInds({{0,1,2,3,4,5,6,7}});
     moris::Matrix< moris::DDSTMat >  tParentNodeRanks(1,8,0);
     moris::Matrix< moris::IndexMat > tEdgeAncestry({{4, 2, 6, 0, 5, 3, 7, 1, 8, 10, 11, 9}});
     moris::Matrix< moris::DDSTMat >  tParentEdgeRanks(1,12,1);
     moris::Matrix< moris::IndexMat > tFaceAncestry({{4, 3, 5, 2, 0, 1}});
     moris::Matrix< moris::DDSTMat >  tParentFaceRanks(1,6,2);
     moris::Matrix< moris::IndexMat > tElementAncestry({{0}});


     // Initialize Template
     Mesh_Modification_Template tRegSubTemplate(tElementAncestry(0,0),
                                                0,
                                                tRegSubNodes,
                                                tParentNodeInds,
                                                tParentNodeRanks,
                                                tEdgeAncestry,
                                                tParentEdgeRanks,
                                                tFaceAncestry,
                                                tParentFaceRanks,
                                                TemplateType::REGULAR_SUBDIVISION_HEX8);

     // Initialize child mesh with template
     Child_Mesh tChildMesh(tRegSubTemplate);

    /*
     * Initialize New Node Ids Indices and auxiliary connectivity
     */
     moris::Matrix< moris::IndexMat > tAddedNodeIndices({{18, 15, 21, 16, 19, 17, 20}});
     moris::Matrix< moris::IndexMat > tInterEdgeIndices({{21, 18, 48, 19, 28, 20, 29}});

    tChildMesh.add_node_indices(tAddedNodeIndices);
    tChildMesh.add_node_ids(tFullNodeIds);

    // Mark intersected edges and the corresponding cm node index (in this case corresponds to proc index)
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(1),tInterEdgeIndices(1),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(3),tInterEdgeIndices(3),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(5),tInterEdgeIndices(5),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(0),tInterEdgeIndices(0),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(4),tInterEdgeIndices(4),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(6),tInterEdgeIndices(6),true);
    tChildMesh.add_entity_to_intersect_connectivity(tAddedNodeIndices(2),tInterEdgeIndices(2),true);

    // Perform mesh modification
    tChildMesh.modify_child_mesh(TemplateType::HIERARCHY_TET4);


    /*
     * Compute child element volumes
     */
    moris::Matrix< moris::IndexMat > const & tElemNode = tChildMesh.get_element_to_node();
    real tTotalChildVol = compute_volume_for_multiple_tets(tNodeCoords,tElemNode);
    CHECK(tTotalChildVol==Approx(1.0));
}
}
