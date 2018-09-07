/*
 * cl_XTK_Cut_Mesh_Reg_to_Nh.cpp
 *
 *  Created on: Sep 25, 2017
 *      Author: ktdoble
 */

#include "catch.hpp"

#include "xtk/cl_XTK_Cut_Mesh.hpp"

// XTKL: Linear Algebra Includes

#include "linalg/cl_XTK_Matrix_Base.hpp"


// DEBUGGING UTILITY INCLUDES
#include "geomeng/fn_Triangle_Geometry.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Builder_Stk.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Tools.hpp"
#include "mesh/cl_Mesh_Bucket.hpp"
#include "mesh/cl_Mesh_Side_Set_Input.hpp"


namespace xtk
{


TEST_CASE("Child Mesh","[Child_Mesh]")
{
    size_t tMax = std::numeric_limits<size_t>::max();
    size_t tParentElementIndex = 1;

    // Processor local node indices
    moris::Matrix<size_t,Default_Matrix_Integer> tNodeInds({{1,2,3,4,5,6}});

    // Locally indexed
    moris::Matrix<size_t,Default_Matrix_Integer> tElementToNode({{1,2,3,4},
                                                       {2,5,3,4},
                                                       {1,2,3,6},
                                                       {1,4,2,5}});

    // Processor local index
   moris::Matrix<size_t,Default_Matrix_Integer> tElementEdgeParentInds({{0, 1, 2, 3, 4, 5},
                                                              {6, 7, 1, 4, 8, 5},
                                                              {0, 1, 2, 9,10,11},
                                                              {6, 7, 1,10,13,14}});

   moris::Matrix<size_t,Default_Matrix_Integer> tElementEdgeParentRanks(4,6,1); // All rank 1

   // Processor local index
   moris::Matrix<size_t,Default_Matrix_Integer> tElementFaceParentInds({{0, 1, 2, 3},
                                                              {4, 5, 6, 1},
                                                              {0, 7, 8, 9},
                                                              {4,10,11, 7}});

   moris::Matrix<size_t,Default_Matrix_Integer> tElementFaceParentRanks(4,4,2); // All rank 2
   moris::Matrix<size_t,Default_Matrix_Integer> tElementInferfaceSides(4,1,tMax);



   Child_Mesh_Test<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tChildMesh(tParentElementIndex,
                                                                                      tNodeInds,
                                                                                      tElementToNode,
                                                                                      tElementEdgeParentInds,
                                                                                      tElementEdgeParentRanks,
                                                                                      tElementFaceParentInds,
                                                                                      tElementFaceParentRanks,
                                                                                      tElementInferfaceSides);


   // BEGIN TEMPLATE ---------------------------------------------------------------------------------------
   // Insert a template
   Mesh_Modification_Template<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tMeshTemplate;
   tMeshTemplate.mNumNewElem = 2;
   tMeshTemplate.mNumElemToReplace = 1;
   tMeshTemplate.mElemIndToReplace = 2;
   tMeshTemplate.mParentEdgeInds           = moris::Matrix<size_t,Default_Matrix_Integer>({{2,4,5,7,8,9},
                                                                                 {0,2,1,4,3,9}});
   tMeshTemplate.mParentEdgeRanks          = moris::Matrix<size_t,Default_Matrix_Integer>(2,6,1);
   tMeshTemplate.mParentFaceInds           = moris::Matrix<size_t,Default_Matrix_Integer>({{2,4,5,7},
                                                                                 {0,2,1,4}});
   tMeshTemplate.mParentFaceRanks          = moris::Matrix<size_t,Default_Matrix_Integer>(2,4,2);
   tMeshTemplate.mNodeInds                 = moris::Matrix<size_t,Default_Matrix_Integer>({{10,11,12,13,14,15,16}});

   tMeshTemplate.mNewElementToNode         = moris::Matrix<size_t,Default_Matrix_Integer>({{0,1,2,3},
                                                                                 {4,5,6,1}});
   tMeshTemplate.mNewParentEdgeRanks       = moris::Matrix<size_t,Default_Matrix_Integer>(2,6,1);

   tMeshTemplate.mNewParentEdgeOrdinals    = moris::Matrix<size_t,Default_Matrix_Integer>({{0,1,2,3,4,5},
                                                                                 {0,1,2,3,4,5}});
   tMeshTemplate.mNewParentFaceRanks       = moris::Matrix<size_t,Default_Matrix_Integer>(2,4,2);

   tMeshTemplate.mNewParentFaceOrdinals    = moris::Matrix<size_t,Default_Matrix_Integer>({{0,1,2,3},
                                                                                 {3,2,1,0}});
   tMeshTemplate.mNewElementInterfaceSides = moris::Matrix<size_t,Default_Matrix_Integer>({{1},
                                                                             {tMax}});
      // END TEMPLATE -----------------------------------------------------------------------------------------

   tChildMesh.add_node_indices(tMeshTemplate.mNodeInds);

   // Insert template into mesh
   tChildMesh.allocate_more_elements(tMeshTemplate.mNumNewElem - tMeshTemplate.mNumElemToReplace);
   tChildMesh.insert_child_mesh_template(tMeshTemplate);
   tChildMesh.generate_connectivities(true,true,true);
}


/*
 * Test the transition between a regular subdivision hex to conformal node hierarchy mesh independently of the geometry engine
 */
TEST_CASE("Regular Subdivision to NH Transition","[RS_NH]")
{
    moris::Matrix<size_t,Default_Matrix_Integer> tRegSubNodes({{0 ,4 ,6 ,2 ,1 ,5 ,7 ,3 ,8 ,9 ,10,11,12,13,14}});
    moris::Matrix<size_t,Default_Matrix_Integer> tRegSubNodeIds({{1,2,4,3,5,6,8,7,9,10,11,12,13,14,15}});
    moris::Matrix<size_t,Default_Matrix_Integer> tFullNodeIds({{1,2 ,4 ,3 ,5 ,6 ,8 ,7 ,9 ,10,11,12,13,14,15,19,16,22,17,20,18,21}});

    /*
     * Setup node coordinate list corresponding to the edge indices above
     */
    moris::Matrix<real,Default_Matrix_Real> tNodeCoords({{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}, {0.5, 0, 0.5}, {1, 0.5, 0.5}, {0.5, 1, 0.5}, {0, 0.5, 0.5}, {0.5, 0.5, 0}, {0.5, 0.5, 1}, {0.5, 0.5, 0.5}, {1, 1, 0.6975}, {1, 0.6975, 1}, {0.6975, 1, 1}, {1, 0.6975, 0.6975}, {0.6975, 1, 0.6975}, {0.6975, 0.6975, 1}, {0.7983333333333332, 0.7983333333333332, 0.7983333333333332}});


    /*
     * Inheritance
     */
     moris::Matrix<size_t,Default_Matrix_Integer> tNodeAncestry({{0}});
     moris::Matrix<size_t,Default_Matrix_Integer> tEdgeAncestry({{4, 2, 6, 0, 5, 3, 7, 1, 8, 10, 11, 9}});
     moris::Matrix<size_t,Default_Matrix_Integer> tParentEdgeRanks(1,12,1);
     moris::Matrix<size_t,Default_Matrix_Integer> tFaceAncestry({{4, 3, 5, 2, 0, 1}});
     moris::Matrix<size_t,Default_Matrix_Integer> tParentFaceRanks(1,6,2);
     moris::Matrix<size_t,Default_Matrix_Integer> tElementAncestry({{0}});


     // Initialize Template
     Mesh_Modification_Template<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tRegSubTemplate(tElementAncestry(0,0),
                                                                                                        0,
                                                                                                        tRegSubNodes,
                                                                                                        tEdgeAncestry,
                                                                                                        tParentEdgeRanks,
                                                                                                        tFaceAncestry,
                                                                                                        tParentFaceRanks,
                                                                                                        TemplateType::REGULAR_SUBDIVISION_HEX8);

     // Initialize child mesh with template
     Child_Mesh_Test<real,size_t,Default_Matrix_Real,Default_Matrix_Integer> tChildMesh(tRegSubTemplate);


    /*
     * Initialize New Node Ids Indices and auxiliary connectivity
     */
     moris::Matrix<size_t,Default_Matrix_Integer> tAddedNodeIndices({{18,
        15,
        21,
        16,
        19,
        17,
        20}});
    moris::Matrix<size_t,Default_Matrix_Integer> tAuxConn(
            {{0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {3, 15, 16, 17,  0,  0,  0, 18, 19, 20,  0,  0,  0},
             {3, 15, 17, 18,  0,  0,  0, 18, 20, 21,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {3, 16, 17, 19,  0,  0,  0, 19, 20, 28,  0,  0,  0},
             {3, 17, 19, 20,  0,  0,  0, 20, 28, 29,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
             {3, 17, 18, 21,  0,  0,  0, 20, 21, 48,  0,  0,  0},
             {3, 17, 20, 21,  0,  0,  0, 20, 29, 48,  0,  0,  0},
             {0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0}});


    tChildMesh.add_node_indices(tAddedNodeIndices);
    tChildMesh.add_node_ids(tFullNodeIds);
    tChildMesh.set_intersect_connectivity(tAuxConn);

    tChildMesh.modify_child_mesh(TemplateType::HIERARCHY_TET4);


    /*
     * Compute child element volumes
     */
    moris::Matrix<size_t,Default_Matrix_Integer> const & tElemNode = tChildMesh.get_element_to_node();
    real tTotalChildVol = compute_volume_for_multiple_tets(tNodeCoords,tElemNode);
    CHECK(tTotalChildVol==Approx(1.0));
}
}
