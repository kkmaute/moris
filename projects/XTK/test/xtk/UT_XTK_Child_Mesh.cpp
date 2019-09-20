/*
 * cl_XTK_Cut_Mesh_Reg_to_Nh.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "cl_XTK_Cut_Mesh.hpp"
#include "xtk/cl_XTK_Child_Mesh_Modification_Template.hpp"

// XTKL: Linear Algebra Includes

#include "cl_Matrix.hpp"


// DEBUGGING UTILITY INCLUDES
#include "geomeng/fn_Triangle_Geometry.hpp"


namespace xtk
{


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

     tChildMesh.add_new_geometry_interface(0);

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
