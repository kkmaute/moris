/*
 * fn_create_edges_from_element_to_node.hpp
 *
 *  Created on: Jun 21, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_
#define SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_


#include "assert/fn_xtk_assert.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"

namespace xtk
{

/*
 * From an element to node connectivity, generate face to node, node to face, element to face connectivity.
 *
 */
template<typename Integer>
void
create_edges_from_element_to_node(enum EntityTopology                 aElementTopology,
                                  Integer                             aNumNodes,
                                  moris::Matrix< moris::IdMat > const & aElementToNode,
                                  moris::Matrix< moris::IdMat >       & aElementToEdge,
                                  moris::Matrix< moris::IdMat >       & aEdgeToNode,
                                  moris::Matrix< moris::IdMat >       & aNodeToEdge,
                                  moris::Matrix< moris::IdMat >       & aEdgeToElement)
{
    XTK_ASSERT(aElementTopology == EntityTopology::TET_4,"This function has only been implemented for tet4 topology");

    //hardcoded values could be provided as a function input
    Integer tMaxEdgePerNode = 100;
    Integer tMaxEdgeToElement = 50;

    // Initialize
    Integer tNumElements = aElementToNode.n_rows();
    Integer tNumEdgePerElem    = 6;
    Integer tNumNodesPerEdge   = 2;
    Integer tNumEdgeCreated    = 0;
    Integer tMaxNumEdges       = tNumElements*tNumEdgePerElem;
    moris::Matrix< moris::IdMat > tNodeToEdgeCounter(1,aNumNodes,0);
    moris::Matrix< moris::IdMat > tEdgeToElemCounter(1,tMaxNumEdges,0);
    Integer tCount = 0;
    Integer tEdgeIndex = 0;
    Integer tNodeInd = 0;
    Integer tFirstInd = 0;

    // Allocate outputs
    aElementToEdge = moris::Matrix< moris::IdMat >(tNumElements,tNumEdgePerElem);
    aEdgeToNode.resize(tMaxNumEdges,tNumNodesPerEdge);
    aNodeToEdge.resize(aNumNodes, tMaxEdgePerNode);
    aNodeToEdge.fill(std::numeric_limits<moris::moris_index>::max());
    aEdgeToElement.resize(tMaxNumEdges,tMaxEdgeToElement);
    aEdgeToElement.fill(std::numeric_limits<moris::moris_index>::max());

    // TET4 specific topology map
    moris::Matrix< moris::IdMat > tElementEdgeToNodeMap({
    {0,1},
    {1,2},
    {0,2},
    {0,3},
    {1,3},
    {2,3}});

    // Single Element Face To Nodes
    moris::Matrix< moris::IdMat > tElementEdgeToNode(6,2);


    Cell<Integer> tPotentialEdges;
    tPotentialEdges.reserve(10);
    Cell<Integer> tPotentialEdges1;
    tPotentialEdges1.reserve(10);
    Cell<Integer> tPotentialEdges2;
    tPotentialEdges2.reserve(10);
    // iterate over elements
    for( Integer i = 0; i<tNumElements; i++)
    {
        tElementEdgeToNode = reindex_matrix(tElementEdgeToNodeMap,i, aElementToNode);

        // iterate over faces in element
        for( Integer j = 0; j<tNumEdgePerElem; j++)
        {

            tEdgeIndex = 0;

            // potential faces (all the ones attached to first node here
            tNodeInd = tElementEdgeToNode(j,tFirstInd);

            // Assemble potential face vector
            Integer tNumPotentialFaces1 = tNodeToEdgeCounter(0,tNodeInd);
            for(Integer k = 0; k< tNumPotentialFaces1; k++)
            {
                tPotentialEdges1.push_back(aNodeToEdge(tNodeInd,k));
            }

            // iterate over nodes on the face j
            for(Integer k = 1; k<tNumNodesPerEdge; k++)
            {
                tNodeInd = tElementEdgeToNode(j,k);

                Integer tNumPotentialFaces2 = tNodeToEdgeCounter(0,tNodeInd);
                for(Integer l = 0; l< tNumPotentialFaces2; l++)
                {
                    tPotentialEdges2.push_back(aNodeToEdge(tNodeInd,l));
                }

                std::set_intersection(tPotentialEdges1.begin(),
                                      tPotentialEdges1.end(),
                                      tPotentialEdges2.begin(),
                                      tPotentialEdges2.end(),
                                      std::back_inserter(tPotentialEdges.data()));

                tPotentialEdges1 = std::move(tPotentialEdges.data());
                tPotentialEdges.clear();
                tPotentialEdges2.clear();

            }

            // If there are no potential faces then create the face
            if(tPotentialEdges1.size() == 0)
            {
                // Add node to face
                for(Integer k = 0; k<tNumNodesPerEdge; k++)
                {
                    tNodeInd = tElementEdgeToNode(j,k);
                    tCount   =  tNodeToEdgeCounter(0,tNodeInd);
                    aNodeToEdge(tNodeInd,tCount) = tNumEdgeCreated;
                    tEdgeIndex = tNumEdgeCreated;
                    tNodeToEdgeCounter(0,tNodeInd)++;
                }
                replace_row(j,tElementEdgeToNode,tNumEdgeCreated,aEdgeToNode);

                tNumEdgeCreated++;
            }

            // if there are two potential faces at this stage that is an issue
            else if(tPotentialEdges1.size() >1)
            {
                std::cout<<"Invalid number of faces found"<<std::endl;
            }
            else
            {
                tEdgeIndex = tPotentialEdges1(0);
            }


            tPotentialEdges1.clear();
            aElementToEdge(i,j) = tEdgeIndex;
            aEdgeToElement(tEdgeIndex,tEdgeToElemCounter(0,tEdgeIndex)) = i;
            tEdgeToElemCounter(0,tEdgeIndex) ++;
        }
    }

    //Remove excess space from output
    aEdgeToNode.resize(tNumEdgeCreated,tNumNodesPerEdge);
    aEdgeToElement.resize(tNumEdgeCreated,tEdgeToElemCounter.max());
}

}



#endif /* SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_ */
