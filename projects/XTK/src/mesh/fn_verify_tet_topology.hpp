/*
 * fn_verify_tet_topology.hpp
 *
 *  Created on: Apr 24, 2018
 *      Author: ktdoble
 */

#ifndef SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_
#define SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_


#include "assert/fn_xtk_assert.hpp"
#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"

namespace xtk
{
template<typename Integer, typename Integer_Matrix>
moris::Mat_New<Integer, Integer_Matrix>
construct_expected_edge_to_node_tet4(Integer iElem,
                                moris::Mat_New<Integer, Integer_Matrix> const & aElementToNode)
{
    // Initialize output
    moris::Mat_New<Integer, Integer_Matrix> tExpectedEdgeToNode(
            {{aElementToNode(iElem,0), aElementToNode(iElem,1)},
             {aElementToNode(iElem,1), aElementToNode(iElem,2)},
             {aElementToNode(iElem,0), aElementToNode(iElem,2)},
             {aElementToNode(iElem,0), aElementToNode(iElem,3)},
             {aElementToNode(iElem,1), aElementToNode(iElem,3)},
             {aElementToNode(iElem,2), aElementToNode(iElem,3)},
            });

    return tExpectedEdgeToNode;
}

template<typename Integer, typename Integer_Matrix>
moris::Mat_New<Integer, Integer_Matrix>
construct_expected_face_to_node_tet4(Integer iElem,
                                moris::Mat_New<Integer, Integer_Matrix> const & aElementToNode)
{
    // construct expected face to node
    moris::Mat_New<Integer, Integer_Matrix> tExpectedFaceToNode(
            {{aElementToNode(iElem,0), aElementToNode(iElem,1), aElementToNode(iElem,3)},
             {aElementToNode(iElem,2), aElementToNode(iElem,1), aElementToNode(iElem,3)},
             {aElementToNode(iElem,0), aElementToNode(iElem,2), aElementToNode(iElem,3)},
             {aElementToNode(iElem,0), aElementToNode(iElem,2), aElementToNode(iElem,1)},
            });

    return tExpectedFaceToNode;
}

template<typename Integer, typename Integer_Matrix>
bool
verify_tet4_edge_topology(moris::Mat_New<Integer, Integer_Matrix> const & aElementToNode,
                          moris::Mat_New<Integer, Integer_Matrix> const & aElementToEdge,
                          moris::Mat_New<Integer, Integer_Matrix> const & aEdgeToNode)
{
    // Number of elements
    bool tValidEdgeTopo = true;

    Integer tNumElem = aElementToNode.n_rows();
    Integer tNumEdgePerElem = aElementToEdge.n_cols();

    // Initilize expected edge to node connectivity
    moris::Mat_New<Integer, Integer_Matrix> tExpectedEdgeToNode;
    moris::Mat_New<Integer, Integer_Matrix> tActualEdgeToNode(6,2);
    moris::Mat_New<Integer, Integer_Matrix> tReorderEdgesMatrix(tNumElem,6);


    XTK_ASSERT(tNumEdgePerElem==6,"TET4 NEEDS TO HAVE 6 EDGES (6 COLUMMNS)");

    for( Integer i = 0; i<tNumElem; i++)
    {

        tExpectedEdgeToNode = construct_expected_edge_to_node_tet4(i,aElementToNode);

        for(Integer j = 0; j<tNumEdgePerElem; j++)
        {

            Integer tEdgeIndex = aElementToEdge(i,j);

            // Make sure the same node does not appear on this edge twice.
            XTK_ASSERT(aEdgeToNode(tEdgeIndex,0)!=aEdgeToNode(tEdgeIndex,1),"Same node cannot appear on the same edge");

            // Loop over nodes in an edge and check for equivalence
            for(Integer k = 0; k<2; k++)
            {

                // Add to actual edge to node
                tActualEdgeToNode(j,k) = aEdgeToNode(tEdgeIndex,k);

                // Check both nodes on edge ( assuming node order on an edge does not matter )
                // Say this is a valid node
                if(   tExpectedEdgeToNode(j,k) ==  aEdgeToNode(tEdgeIndex,0)
                   || tExpectedEdgeToNode(j,k) ==  aEdgeToNode(tEdgeIndex,1) )
                {
                    // Intentionally empty
                }

                // Invalid topo
                else
                {
                    std::cout<<"iElem = "<< i<< " iEdge = "<< j<<std::endl;
                    std::cout<< "Expected Nodes = "<< tExpectedEdgeToNode(j,0)<<" "<<tExpectedEdgeToNode(j,1)<<std::endl;
                    std::cout<< "Actual Nodes   = "<< aEdgeToNode(tEdgeIndex,0)<<" "<<aEdgeToNode(tEdgeIndex,1)<<std::endl;
                    tValidEdgeTopo = false;
                }
            }

        }


        // Reorder edges as follows
        for(Integer iA = 0; iA<tNumEdgePerElem; iA++)
        {
            for (Integer iExp = 0; iExp<tNumEdgePerElem; iExp++)
            {

                bool tN1Found = tExpectedEdgeToNode(iExp,0) ==  tActualEdgeToNode(iA,0) || tExpectedEdgeToNode(iExp,0) ==  tActualEdgeToNode(iA,1);
                bool tN2Found = tExpectedEdgeToNode(iExp,1) ==  tActualEdgeToNode(iA,0) || tExpectedEdgeToNode(iExp,1) ==  tActualEdgeToNode(iA,1);

                if(tN1Found && tN2Found)
                {
                    tReorderEdgesMatrix(i,iA) = iExp;
                }
            }
        }


    }

    if(!tValidEdgeTopo)
    {
    print(tReorderEdgesMatrix,"tReorderEdgesMatrix");
    }
    return tValidEdgeTopo;
}


template<typename Integer, typename Integer_Matrix>
bool
verify_tet4_face_topology(moris::Mat_New<Integer, Integer_Matrix> const & aElementToNode,
                          moris::Mat_New<Integer, Integer_Matrix> const & aElementToFace,
                          moris::Mat_New<Integer, Integer_Matrix> const & aFaceToNode)
{
    // Number of elements
        bool tValidFaceTopo = true;

        Integer tNumElem        = aElementToNode.n_rows();
        Integer tNumFacePerElem = aElementToFace.n_cols();
        Integer tNumNodePerFace = 3;

        // Initilize expected edge to node connectivity
        moris::Mat_New<Integer, Integer_Matrix> tExpectedEdgeToNode;
        moris::Mat_New<Integer, Integer_Matrix> tActualEdgeToNode(6,3);
        moris::Mat_New<Integer, Integer_Matrix> tReorderEdgesMatrix(tNumElem,6);


        XTK_ASSERT(tNumFacePerElem==4,"TET4 NEEDS TO HAVE 4 FACES (4 COLUMMNS)");

        for( Integer i = 0; i<tNumElem; i++)
        {

            tExpectedEdgeToNode = construct_expected_face_to_node_tet4(i,aElementToNode);

            for(Integer j = 0; j<tNumFacePerElem; j++)
            {

                Integer tFaceIndex = aElementToFace(i,j);

                // Make sure the same node does not appear on this face twice.
                XTK_ASSERT(aFaceToNode(tFaceIndex,0)!=aFaceToNode(tFaceIndex,1),"Same node cannot appear on the same edge");
                XTK_ASSERT(aFaceToNode(tFaceIndex,0)!=aFaceToNode(tFaceIndex,2),"Same node cannot appear on the same edge");
                XTK_ASSERT(aFaceToNode(tFaceIndex,1)!=aFaceToNode(tFaceIndex,2),"Same node cannot appear on the same edge");

                // Loop over nodes in an edge and check for equivalence
                for(Integer k = 0; k<tNumNodePerFace; k++)
                {

                    // Add to actual edge to node
                    tActualEdgeToNode(j,k) = aFaceToNode(tFaceIndex,k);

                    // Check all nodes on face ( assuming node order on an edge does not matter )
                    if(   tExpectedEdgeToNode(j,k) ==  aFaceToNode(tFaceIndex,0)
                       || tExpectedEdgeToNode(j,k) ==  aFaceToNode(tFaceIndex,1)
                       || tExpectedEdgeToNode(j,k) ==  aFaceToNode(tFaceIndex,2))
                    {
                        // Intentionally empty
                    }

                    // Invalid topo
                    else
                    {
                        std::cout<<"iElem = "<< i<< " iFace = "<< j<<std::endl;
                        std::cout<< "Expected Nodes = "<< tExpectedEdgeToNode(j,0)<<" "<<tExpectedEdgeToNode(j,1)<<" "<<tExpectedEdgeToNode(j,2)<<std::endl;
                        std::cout<< "Actual Nodes   = "<< aFaceToNode(tFaceIndex,0)  <<" "<<aFaceToNode(tFaceIndex,1)  <<" "<<aFaceToNode(tFaceIndex,2)<<std::endl;
                        tValidFaceTopo = false;
                    }
                }

            }

        }

        return tValidFaceTopo;
}


/*
 * Verify that the tet4 topology provided results in a valid topology. (i.e. edge 0 contains node 0 and node 1)
 * This function checks the face ordering, edge ordering
 * @param[in] aElementNodes -
 */
template <typename Integer, typename Integer_Matrix>
bool
verify_tet4_topology(moris::Mat_New<Integer, Integer_Matrix> const & aElementToNode,
                     moris::Mat_New<Integer, Integer_Matrix> const & aElementToEdge,
                     moris::Mat_New<Integer, Integer_Matrix> const & aElementToFace,
                     moris::Mat_New<Integer, Integer_Matrix> const & aEdgeToNode,
                     moris::Mat_New<Integer, Integer_Matrix> const & aFaceToNode)
{
    XTK_ASSERT(aElementToNode.n_cols() == 4,"INVALID NUMBER OF NODES PROVIDED, NEEDS TO BE 4 FOR TET4");

    bool tValidTopo = false;


    // verify the edges
    bool tValidEdgeTopo = verify_tet4_edge_topology(aElementToNode,aElementToEdge,aEdgeToNode);

    // Verify faces
    bool tValidFaceTopo = verify_tet4_face_topology(aElementToNode,aElementToFace,aFaceToNode);

    if( !tValidEdgeTopo)
    {
        std::cout<<"Invalid edge topology detected"<<std::endl;
        tValidTopo = false;
    }

    else if(!tValidFaceTopo)
    {
        std::cout<<"Invalid face topology detected"<<std::endl;
        tValidTopo = false;
    }

    else
    {
        tValidTopo = true;
    }
    return tValidTopo;

}


}


#endif /* SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_ */
