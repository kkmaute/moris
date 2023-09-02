/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_create_faces_from_element_to_node.hpp
 *
 */

#ifndef SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_
#define SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_

#include <unordered_map>

#include "cl_Mesh_Enums.hpp"
#include "cl_MTK_Cell_Info_Tet4.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "cl_Matrix.hpp"
#include "cl_MTK_Cell_Info.hpp"
namespace xtk
{

    /*
     * From an element to node connectivity, generate face to node, node to face, element to face connectivity.
     *
     */
    inline
    void
    create_faces_from_element_to_node(
            moris::mtk::Cell_Info*                   aCellConnectivity,
            moris::size_t                            aNumNodes,
            moris::Matrix< moris::IndexMat > const & aElementToNode,
            moris::Matrix< moris::IndexMat >       & aElementToFace,
            moris::Matrix< moris::IndexMat >       & aFaceToNode,
            moris::Matrix< moris::IndexMat >       & aNodeToFace,
            moris::Matrix< moris::IndexMat >       & aFaceToElement)
    {
        //hard-coded values could be provided as a function input
        //FIXME: convert aNodeToFace into cell of matrices
        moris::size_t tMaxFacePerNode = 10;
        moris::size_t tMaxUsed        = 0;

        // Initialize
        moris::size_t tNumElements     = aElementToNode.n_rows();
        moris::size_t tNumFacesPerElem = aCellConnectivity->get_num_facets();
        moris::size_t tNumNodesPerFace = aCellConnectivity->get_num_verts_per_facet();
        moris::size_t tNumFaceCreated  = 0;
        moris::size_t tMaxNumFaces     = tNumElements*tNumFacesPerElem;

        moris::Matrix< moris::IndexMat > tNodeToFaceCounter(1,aNumNodes,0);
        moris::Matrix< moris::IndexMat > tFaceToElemCounter(1,tMaxNumFaces,0);

        moris::size_t tCount     = 0;
        moris::size_t tFaceIndex = 0;
        moris::size_t tNodeInd   = 0;
        moris::size_t tFirstInd  = 0;
        moris::size_t tResize    = 0;

        // Allocate outputs
        aElementToFace.resize(tNumElements,tNumFacesPerElem);
        aFaceToNode.resize   (tMaxNumFaces,tNumNodesPerFace);

        aNodeToFace.resize(aNumNodes, tMaxFacePerNode);
        aNodeToFace.fill(MORIS_INDEX_MAX);

        aFaceToElement.resize(tMaxNumFaces,2);
        aFaceToElement.fill(MORIS_INDEX_MAX);

        // TET4 specific topology map
        moris::Matrix< moris::IndexMat > tNodeToFaceMap =  aCellConnectivity->get_node_to_facet_map();

        // Single Element Face To Nodes
        moris::Matrix< moris::IndexMat > tElementFaceToNode(1,tNumNodesPerFace);

        moris::Cell<moris::size_t> tPotentialFaces;
        moris::Cell<moris::size_t> tPotentialFaces1;
        moris::Cell<moris::size_t> tPotentialFaces2;

        // iterate over elements
        for( moris::size_t i = 0; i<tNumElements; i++)
        {
            tElementFaceToNode = reindex_matrix(tNodeToFaceMap,i, aElementToNode);

            // iterate over faces in element
            for( moris::size_t j = 0; j<tNumFacesPerElem; j++)
            {
                tFaceIndex = 0;

                // potential faces (all the ones attached to first node here
                tNodeInd = tElementFaceToNode(j,tFirstInd);

                // Assemble potential face vector
                moris::size_t tNumPotentialFaces1 = tNodeToFaceCounter(0,tNodeInd);
                for(moris::size_t k = 0; k< tNumPotentialFaces1; k++)
                {
                    tPotentialFaces1.push_back(aNodeToFace(tNodeInd,k));
                }

                // iterate over nodes on the face j
                for(moris::size_t k = 1; k<tNumNodesPerFace; k++)
                {
                    tNodeInd = tElementFaceToNode(j,k);

                    moris::size_t tNumPotentialFaces2 = tNodeToFaceCounter(0,tNodeInd);
                    for(moris::size_t l = 0; l< tNumPotentialFaces2; l++)
                    {
                        tPotentialFaces2.push_back(aNodeToFace(tNodeInd,l));
                    }

                    std::set_intersection(tPotentialFaces1.begin(),
                            tPotentialFaces1.end(),
                            tPotentialFaces2.begin(),
                            tPotentialFaces2.end(),
                            std::back_inserter(tPotentialFaces));

                    tPotentialFaces1 = std::move(tPotentialFaces.data());
                    tPotentialFaces.clear();
                    tPotentialFaces2.clear();
                }

                // If there are no potential faces then create the face
                if(tPotentialFaces1.size() == 0)
                {
                    // Add node to face
                    for(moris::size_t k = 0; k<tNumNodesPerFace; k++)
                    {
                        tNodeInd = tElementFaceToNode(j,k);
                        tCount   = tNodeToFaceCounter(0,tNodeInd);

                        // make sure we have not exceeded the allocated space in node to face
                        if(tCount>=aNodeToFace.n_cols())
                        {
                            aNodeToFace.resize(aNumNodes,aNodeToFace.n_cols()+tMaxFacePerNode);
                            tResize++;
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

        // Check for number of resize operations; note check needs to be "<=" for tNumElements = tResize = 0
        MORIS_CHECK_MEMORY( tResize <= tNumElements,
                "create_faces_from_element_to_node: Number of resize operations too large (%zu / %zu) - increase tMaxFacePerNode parameter.\n",
                tResize,tNumElements);
    }
}

#endif /* SRC_XTK_FN_CREATE_FACES_FROM_ELEMENT_TO_NODE_HPP_ */

