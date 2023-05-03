/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_create_edges_from_element_to_node.hpp
 *
 */

#ifndef SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_
#define SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_

#include "cl_Mesh_Enums.hpp"
#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"

namespace xtk
{

    /*
     * From an element to node connectivity, generate face to node, node to face, element to face connectivity.
     *
     */
    template<typename Integer>
    inline void
    create_edges_from_element_to_node(
            enum CellTopology                   aElementTopology,
            Integer                             aNumNodes,
            moris::Matrix< moris::IdMat > const & aElementToNode,
            moris::Matrix< moris::IdMat >       & aElementToEdge,
            moris::Matrix< moris::IdMat >       & aEdgeToNode,
            moris::Matrix< moris::IdMat >       & aNodeToEdge,
            moris::Matrix< moris::IdMat >       & aEdgeToElement)
    {
        //hard-coded values could be provided as a function input
        Integer tMaxEdgePerNode   = 20;
        Integer tMaxEdgeToElement = 10;

        // Initialize
        Integer tNumEdgePerElem = 0;

        moris::Matrix< moris::IdMat > tElementEdgeToNodeMap(0,0);

        switch(aElementTopology)
        {
            case CellTopology::TET4:
            {
                tNumEdgePerElem = 6;

                // TET4 specific topology map
                tElementEdgeToNodeMap = {
                        {0,1},
                        {1,2},
                        {0,2},
                        {0,3},
                        {1,3},
                        {2,3}};

                break;
            }
            case CellTopology::TRI3:
            {
                tNumEdgePerElem = 3;

                // QUAD4 specific topology map
                tElementEdgeToNodeMap = {
                        {0,1},
                        {1,2},
                        {0,2}};

                break;
            }
            default:
            {
                MORIS_ASSERT(0, "This function has only been implemented for tet4 and tri3 topology.");
            }
        }

        Integer tNumElements       = aElementToNode.n_rows();
        Integer tNumNodesPerEdge   = 2;
        Integer tNumEdgeCreated    = 0;
        Integer tMaxNumEdges       = tNumElements*tNumEdgePerElem;

        moris::Matrix< moris::IdMat > tNodeToEdgeCounter(1,aNumNodes,0);
        moris::Matrix< moris::IdMat > tEdgeToElemCounter(1,tMaxNumEdges,0);

        Integer tCount     = 0;
        Integer tEdgeIndex = 0;
        Integer tNodeInd   = 0;
        Integer tFirstInd  = 0;
        Integer tNumResize = 0;

        // Allocate outputs
        aElementToEdge = moris::Matrix< moris::IdMat >(tNumElements,tNumEdgePerElem);

        aEdgeToNode.resize(tMaxNumEdges,tNumNodesPerEdge);
        aNodeToEdge.resize(aNumNodes, tMaxEdgePerNode);
        aEdgeToElement.resize(tMaxNumEdges,tMaxEdgeToElement);

        aNodeToEdge.fill(std::numeric_limits<moris::moris_index>::max());
        aEdgeToElement.fill(std::numeric_limits<moris::moris_index>::max());

        // Single Element Face To Nodes
        moris::Matrix< moris::IdMat > tElementEdgeToNode(tNumEdgePerElem,tNumNodesPerEdge);

        moris::Cell<Integer> tPotentialEdges;
        tPotentialEdges.reserve(10);
        moris::Cell<Integer> tPotentialEdges1;
        tPotentialEdges1.reserve(10);
        moris::Cell<Integer> tPotentialEdges2;
        tPotentialEdges2.reserve(10);

        // iterate over elements
        for( Integer i = 0; i<tNumElements; i++)
        {
            tElementEdgeToNode = reindex_matrix(tElementEdgeToNodeMap,i, aElementToNode);

            // iterate over edges in element
            for( Integer j = 0; j<tNumEdgePerElem; j++)
            {
                tEdgeIndex = 0;

                // potential edges (all the ones attached to first node here)
                tNodeInd = tElementEdgeToNode(j,tFirstInd);

                // Assemble potential edge vector
                Integer tNumPotentialEdges1 = tNodeToEdgeCounter(0,tNodeInd);
                for(Integer k = 0; k< tNumPotentialEdges1; k++)
                {
                    tPotentialEdges1.push_back(aNodeToEdge(tNodeInd,k));
                }

                // iterate over nodes on the edge j
                for(Integer k = 1; k<tNumNodesPerEdge; k++)
                {
                    tNodeInd = tElementEdgeToNode(j,k);

                    Integer tNumPotentialEdges2 = tNodeToEdgeCounter(0,tNodeInd);
                    for(Integer l = 0; l< tNumPotentialEdges2; l++)
                    {
                        tPotentialEdges2.push_back(aNodeToEdge(tNodeInd,l));
                    }

                    std::set_intersection(tPotentialEdges1.begin(),
                            tPotentialEdges1.end(),
                            tPotentialEdges2.begin(),
                            tPotentialEdges2.end(),
                            std::back_inserter(tPotentialEdges));

                    tPotentialEdges1 = std::move(tPotentialEdges.data());
                    tPotentialEdges.clear();
                    tPotentialEdges2.clear();
                }

                // If there are no potential edges then create the edge
                if(tPotentialEdges1.size() == 0)
                {
                    // Add node to edge
                    for(Integer k = 0; k<tNumNodesPerEdge; k++)
                    {
                        tNodeInd   = tElementEdgeToNode(j,k);
                        tCount     = tNodeToEdgeCounter(0,tNodeInd);
                        tEdgeIndex = tNumEdgeCreated;

                        if (tCount >= aNodeToEdge.n_cols() )
                        {
                            aNodeToEdge.resize( aNumNodes, aNodeToEdge.n_cols() + tMaxEdgePerNode );
                            tNumResize++;
                        }

                        aNodeToEdge(tNodeInd,tCount) = tNumEdgeCreated;

                        tNodeToEdgeCounter(0,tNodeInd)++;
                    }
                    replace_row(j,tElementEdgeToNode,tNumEdgeCreated,aEdgeToNode);

                    tNumEdgeCreated++;
                }

                // if there are two potential edges at this stage that is an issue
                else if(tPotentialEdges1.size() >1)
                {
                    std::cout<<"Invalid number of edges found"<<std::endl;
                }
                else
                {
                    tEdgeIndex = tPotentialEdges1(0);
                }

                tPotentialEdges1.clear();
                aElementToEdge(i,j) = tEdgeIndex;
                tCount = tEdgeToElemCounter(0,tEdgeIndex);

                if ( tCount >= aEdgeToElement.n_cols() )
                {
                    aEdgeToElement.resize( tMaxNumEdges, aEdgeToElement.n_cols() + tMaxEdgeToElement );
                    tNumResize++;
                }

                aEdgeToElement(tEdgeIndex,tCount) = i;
                tEdgeToElemCounter(0,tEdgeIndex) ++;
            }
        }

        MORIS_CHECK_MEMORY(tNumResize < 2*tNumElements,
                "Excessive use of resize; increase parameter tMaxEdgePerNode: %d times for %d elements.\n",
                tNumResize,tNumElements);

        //Remove excess space from output
        aNodeToEdge.resize(aNumNodes,tNodeToEdgeCounter.max());
        aEdgeToNode.resize(tNumEdgeCreated,tNumNodesPerEdge);
        aEdgeToElement.resize(tNumEdgeCreated,tEdgeToElemCounter.max());
    }
}

#endif /* SRC_XTK_FN_CREATE_EDGES_FROM_ELEMENT_TO_NODE_HPP_ */

