/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_verify_tet_topology.hpp
 *
 */

#ifndef SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_
#define SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "fn_print.hpp"

namespace moris::xtk
{
    inline Matrix< IndexMat >
    construct_expected_edge_to_node_tet4(
            moris::moris_index        iElem,
            Matrix< IndexMat > const &aElementToNode )
    {
        // Initialize output
        Matrix< IndexMat > tExpectedEdgeToNode(
                {
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 1 ) },
                        { aElementToNode( iElem, 1 ), aElementToNode( iElem, 2 ) },
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 2 ) },
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 3 ) },
                        { aElementToNode( iElem, 1 ), aElementToNode( iElem, 3 ) },
                        { aElementToNode( iElem, 2 ), aElementToNode( iElem, 3 ) },
                } );

        return tExpectedEdgeToNode;
    }

    //------------------------------------------------------------------------------

    inline Matrix< IndexMat >
    construct_expected_face_to_node_tet4(
            moris::moris_index        iElem,
            Matrix< IndexMat > const &aElementToNode )
    {
        // construct expected face to node
        Matrix< IndexMat > tExpectedFaceToNode(
                {
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 1 ), aElementToNode( iElem, 3 ) },
                        { aElementToNode( iElem, 2 ), aElementToNode( iElem, 1 ), aElementToNode( iElem, 3 ) },
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 2 ), aElementToNode( iElem, 3 ) },
                        { aElementToNode( iElem, 0 ), aElementToNode( iElem, 2 ), aElementToNode( iElem, 1 ) },
                } );

        return tExpectedFaceToNode;
    }

    //------------------------------------------------------------------------------

    inline bool
    verify_tet4_edge_topology(
            Matrix< IndexMat > const &aElementToNode,
            Matrix< IndexMat > const &aElementToEdge,
            Matrix< IndexMat > const &aEdgeToNode )
    {

        // Number of elements
        bool tValidEdgeTopo = true;

        moris::moris_index tNumElem        = aElementToNode.n_rows();
        moris::moris_index tNumEdgePerElem = aElementToEdge.n_cols();

        // Initilize expected edge to node connectivity
        Matrix< IndexMat > tExpectedEdgeToNode;
        Matrix< IndexMat > tActualEdgeToNode( 6, 2 );
        Matrix< IndexMat > tReorderEdgesMatrix( tNumElem, 6 );

        MORIS_ASSERT( tNumEdgePerElem == 6, "TET4 NEEDS TO HAVE 6 EDGES (6 COLUMMNS)" );

        for ( moris::moris_index i = 0; i < tNumElem; i++ )
        {

            tExpectedEdgeToNode = construct_expected_edge_to_node_tet4( i, aElementToNode );

            for ( moris::moris_index j = 0; j < tNumEdgePerElem; j++ )
            {

                moris::moris_index tEdgeIndex = aElementToEdge( i, j );

                // Make sure the same node does not appear on this edge twice.
                MORIS_ASSERT( aEdgeToNode( tEdgeIndex, 0 ) != aEdgeToNode( tEdgeIndex, 1 ), "Same node cannot appear on the same edge" );

                // Loop over nodes in an edge and check for equivalence
                for ( moris::moris_index k = 0; k < 2; k++ )
                {

                    // Add to actual edge to node
                    tActualEdgeToNode( j, k ) = aEdgeToNode( tEdgeIndex, k );

                    // Check both nodes on edge ( assuming node order on an edge does not matter )
                    // Say this is a valid node
                    if ( tExpectedEdgeToNode( j, k ) == aEdgeToNode( tEdgeIndex, 0 )
                            || tExpectedEdgeToNode( j, k ) == aEdgeToNode( tEdgeIndex, 1 ) )
                    {
                        // Intentionally empty
                    }

                    // Invalid topo
                    else
                    {
                        std::cout << "iElem = " << i << " iEdge = " << j << '\n';
                        std::cout << "Expected Nodes = " << tExpectedEdgeToNode( j, 0 ) << " " << tExpectedEdgeToNode( j, 1 ) << '\n';
                        std::cout << "Actual Nodes   = " << aEdgeToNode( tEdgeIndex, 0 ) << " " << aEdgeToNode( tEdgeIndex, 1 ) << '\n';
                        tValidEdgeTopo = false;
                    }
                }
            }

            // Reorder edges as follows
            for ( moris::moris_index iA = 0; iA < tNumEdgePerElem; iA++ )
            {
                for ( moris::moris_index iExp = 0; iExp < tNumEdgePerElem; iExp++ )
                {

                    bool tN1Found = tExpectedEdgeToNode( iExp, 0 ) == tActualEdgeToNode( iA, 0 ) || tExpectedEdgeToNode( iExp, 0 ) == tActualEdgeToNode( iA, 1 );
                    bool tN2Found = tExpectedEdgeToNode( iExp, 1 ) == tActualEdgeToNode( iA, 0 ) || tExpectedEdgeToNode( iExp, 1 ) == tActualEdgeToNode( iA, 1 );

                    if ( tN1Found && tN2Found )
                    {
                        tReorderEdgesMatrix( i, iA ) = iExp;
                    }
                }
            }
        }

        if ( !tValidEdgeTopo )
        {
            moris::print( tReorderEdgesMatrix, "tReorderEdgesMatrix" );
        }
        return tValidEdgeTopo;
    }

    //------------------------------------------------------------------------------

    inline bool
    verify_tet4_face_topology(
            Matrix< IndexMat > const &aElementToNode,
            Matrix< IndexMat > const &aElementToFace,
            Matrix< IndexMat > const &aFaceToNode )
    {
        // Number of elements
        bool tValidFaceTopo = true;

        moris::moris_index tNumElem        = aElementToNode.n_rows();
        moris::moris_index tNumFacePerElem = aElementToFace.n_cols();
        moris::moris_index tNumNodePerFace = 3;

        // Initilize expected edge to node connectivity
        Matrix< IndexMat > tExpectedEdgeToNode;
        Matrix< IndexMat > tActualEdgeToNode( 6, 3 );
        Matrix< IndexMat > tReorderEdgesMatrix( tNumElem, 6 );

        MORIS_ASSERT( tNumFacePerElem == 4, "TET4 NEEDS TO HAVE 4 FACES (4 COLUMMNS)" );

        for ( moris::moris_index i = 0; i < tNumElem; i++ )
        {

            tExpectedEdgeToNode = construct_expected_face_to_node_tet4( i, aElementToNode );

            for ( moris::moris_index j = 0; j < tNumFacePerElem; j++ )
            {

                moris::moris_index tFaceIndex = aElementToFace( i, j );

                // Make sure the same node does not appear on this face twice.
                MORIS_ASSERT( aFaceToNode( tFaceIndex, 0 ) != aFaceToNode( tFaceIndex, 1 ), "Same node cannot appear on the same edge" );
                MORIS_ASSERT( aFaceToNode( tFaceIndex, 0 ) != aFaceToNode( tFaceIndex, 2 ), "Same node cannot appear on the same edge" );
                MORIS_ASSERT( aFaceToNode( tFaceIndex, 1 ) != aFaceToNode( tFaceIndex, 2 ), "Same node cannot appear on the same edge" );

                // Loop over nodes in an edge and check for equivalence
                for ( moris::moris_index k = 0; k < tNumNodePerFace; k++ )
                {

                    // Add to actual edge to node
                    tActualEdgeToNode( j, k ) = aFaceToNode( tFaceIndex, k );

                    // Check all nodes on face ( assuming node order on an edge does not matter )
                    if ( tExpectedEdgeToNode( j, k ) == aFaceToNode( tFaceIndex, 0 )
                            || tExpectedEdgeToNode( j, k ) == aFaceToNode( tFaceIndex, 1 )
                            || tExpectedEdgeToNode( j, k ) == aFaceToNode( tFaceIndex, 2 ) )
                    {
                        // Intentionally empty
                    }

                    // Invalid topo
                    else
                    {
                        std::cout << "iElem = " << i << " iFace = " << j << '\n';
                        std::cout << "Expected Nodes = " << tExpectedEdgeToNode( j, 0 ) << " " << tExpectedEdgeToNode( j, 1 ) << " " << tExpectedEdgeToNode( j, 2 ) << '\n';
                        std::cout << "Actual Nodes   = " << aFaceToNode( tFaceIndex, 0 ) << " " << aFaceToNode( tFaceIndex, 1 ) << " " << aFaceToNode( tFaceIndex, 2 ) << '\n';
                        tValidFaceTopo = false;
                    }
                }
            }
        }

        return tValidFaceTopo;
    }

    //------------------------------------------------------------------------------

    /*
     * Verify that the tet4 topology provided results in a valid topology. (i.e. edge 0 contains node 0 and node 1)
     * This function checks the face ordering, edge ordering
     * @param[in] aElementNodes -
     */

    inline bool
    verify_tet4_topology(
            Matrix< IndexMat > const &aElementToNode,
            Matrix< IndexMat > const &aElementToEdge,
            Matrix< IndexMat > const &aElementToFace,
            Matrix< IndexMat > const &aEdgeToNode,
            Matrix< IndexMat > const &aFaceToNode )
    {
        MORIS_ASSERT( aElementToNode.n_cols() == 4, "UNDEFINED NUMBER OF NODES PROVIDED, NEEDS TO BE 4 FOR TET4" );

        bool tValidTopo = false;

        // verify the edges
        bool tValidEdgeTopo = verify_tet4_edge_topology( aElementToNode, aElementToEdge, aEdgeToNode );

        // Verify faces
        bool tValidFaceTopo = verify_tet4_face_topology( aElementToNode, aElementToFace, aFaceToNode );

        if ( !tValidEdgeTopo )
        {
            std::cout << "Invalid edge topology detected" << '\n';
            tValidTopo = false;
        }

        else if ( !tValidFaceTopo )
        {
            std::cout << "Invalid face topology detected" << '\n';
            tValidTopo = false;
        }

        else
        {
            tValidTopo = true;
        }
        return tValidTopo;
    }

}    // namespace moris::xtk

#endif /* SRC_MESH_FN_VERIFY_TET_TOPOLOGY_HPP_ */
