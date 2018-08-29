/*
 * fn_create_faces.hpp
 *
 *  Created on: Jun 19, 2018
 *      Author: ktdoble
 */

#ifndef SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_
#define SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_

#include <unordered_map>
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
template<typename Integer, typename Integer_Matrix>
void
create_faces_from_element_to_node(enum EntityTopology                         aElementTopology,
                                  Integer                                     aNumNodes,
                                  moris::Matrix<Integer, Integer_Matrix> const & aElementToNode,
                                  moris::Matrix<Integer, Integer_Matrix>       & aElementToFace,
                                  moris::Matrix<Integer, Integer_Matrix>       & aFaceToNode,
                                  moris::Matrix<Integer, Integer_Matrix>       & aNodeToFace,
                                  moris::Matrix<Integer, Integer_Matrix>       & aFaceToElement)
{
    XTK_ASSERT(aElementTopology == EntityTopology::TET_4,"This function has only been tested with tet4 topology");

    //hardcoded values could be provided as a function input
    Integer tMaxFacePerNode = 10;
    Integer tMaxUsed = 0;

    // Initialize
    Integer tNumElements = aElementToNode.n_rows();
    Integer tNumFacesPerElem   = 4;
    Integer tNumNodesPerFace   = 3;
    Integer tNumFaceCreated    = 0;
    Integer tMaxNumFaces       = tNumElements*tNumFacesPerElem;
    moris::Matrix<Integer, Integer_Matrix> tNodeToFaceCounter(1,aNumNodes,0);
    moris::Matrix<Integer, Integer_Matrix> tFaceToElemCounter(1,tMaxNumFaces,0);
    Integer tCount = 0;
    Integer tFaceIndex = 0;
    Integer tNodeInd = 0;
    Integer tFirstInd = 0;

    // Allocate outputs
    aElementToFace.resize(tNumElements,tNumFacesPerElem);
    aFaceToNode.resize(tMaxNumFaces,tNumNodesPerFace);
    aNodeToFace.resize(aNumNodes, tMaxFacePerNode);
    aNodeToFace.fill(std::numeric_limits<Integer>::max());
    aFaceToElement.resize(tMaxNumFaces,10);
    aFaceToElement.fill(std::numeric_limits<Integer>::max());

    // TET4 specific topology map
    moris::Matrix<Integer, Integer_Matrix> tElementFacesToNodeMap(
                            {{0, 1, 3},
                             {2, 1, 3},
                             {0, 2, 3},
                             {0, 2, 1}});

    // Single Element Face To Nodes
    moris::Matrix<Integer, Integer_Matrix> tElementFaceToNode;


    Cell<Integer> tPotentialFaces;
    tPotentialFaces.reserve(10);
    Cell<Integer> tPotentialFaces1;
    tPotentialFaces1.reserve(10);
    Cell<Integer> tPotentialFaces2;
    tPotentialFaces2.reserve(10);
    // iterate over elements
    for( Integer i = 0; i<tNumElements; i++)
    {
        tElementFaceToNode = reindex_matrix(tElementFacesToNodeMap,i, aElementToNode);

        // iterate over faces in element
        for( Integer j = 0; j<tNumFacesPerElem; j++)
        {

            tFaceIndex = 0;

            // potential faces (all the ones attached to first node here
            tNodeInd = tElementFaceToNode(j,tFirstInd);

            // Assemble potential face vector
            Integer tNumPotentialFaces1 = tNodeToFaceCounter(0,tNodeInd);
            for(Integer k = 0; k< tNumPotentialFaces1; k++)
            {
                tPotentialFaces1.push_back(aNodeToFace(tNodeInd,k));
            }

            // iterate over nodes on the face j
            for(Integer k = 1; k<tNumNodesPerFace; k++)
            {
                tNodeInd = tElementFaceToNode(j,k);

                Integer tNumPotentialFaces2 = tNodeToFaceCounter(0,tNodeInd);
                for(Integer l = 0; l< tNumPotentialFaces2; l++)
                {
                    tPotentialFaces2.push_back(aNodeToFace(tNodeInd,l));
                }

                std::set_intersection(tPotentialFaces1.begin(),
                                      tPotentialFaces1.end(),
                                      tPotentialFaces2.begin(),
                                      tPotentialFaces2.end(),
                                      std::back_inserter(tPotentialFaces.data()));

                tPotentialFaces1 = std::move(tPotentialFaces.data());
                tPotentialFaces.clear();
                tPotentialFaces2.clear();

            }

            // If there are no potential faces then create the face
            if(tPotentialFaces1.size() == 0)
            {
                // Add node to face
                for(Integer k = 0; k<tNumNodesPerFace; k++)
                {
                    tNodeInd = tElementFaceToNode(j,k);
                    tCount =  tNodeToFaceCounter(0,tNodeInd);

                    // make sure we havent exceeded the allocatd space in node to face
                    if(tCount>=aNodeToFace.n_cols())
                    {
                        aNodeToFace.resize(aNumNodes,aNodeToFace.n_cols()+tMaxFacePerNode);
                    }

                    if(tCount>tMaxUsed)
                    {
                        tMaxUsed=tCount;
                    }
                    aNodeToFace(tNodeInd,tCount) = tNumFaceCreated;
                    tFaceIndex = tNumFaceCreated;
                    tNodeToFaceCounter(0,tNodeInd)++;
                }
                replace_row(j,tElementFaceToNode,tNumFaceCreated,aFaceToNode);

                tNumFaceCreated++;
            }

            // if there are two potential faces at this stage that is an issue
            else if(tPotentialFaces1.size() >1)
            {
                std::cout<<"Invalid number of faces found"<<std::endl;
            }
            else
            {
                tFaceIndex = tPotentialFaces1(0);
            }


            tPotentialFaces1.clear();
            aElementToFace(i,j) = tFaceIndex;
            aFaceToElement(tFaceIndex,tFaceToElemCounter(0,tFaceIndex)) = i;
            tFaceToElemCounter(0,tFaceIndex) ++;
        }
    }

    //Remove excess space from output
    aFaceToNode.resize(tNumFaceCreated,tNumNodesPerFace);
    aFaceToElement.resize(tNumFaceCreated,2);
    aNodeToFace.resize(aNumNodes,tMaxUsed);
}




}



#endif /* SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_ */
