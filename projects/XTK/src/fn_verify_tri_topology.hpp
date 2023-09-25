/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_verify_tri_topology.hpp
 *
 */

#ifndef SRC_MESH_FN_VERIFY_TRI_TOPOLOGY_HPP_
#define SRC_MESH_FN_VERIFY_TRI_TOPOLOGY_HPP_

#include "cl_Matrix.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "fn_print.hpp"

namespace xtk
{

    inline moris::Matrix< moris::IndexMat >
    construct_expected_edge_to_node_tri3(
            moris::moris_index                       iElem,
            moris::Matrix< moris::IndexMat > const & aElementToNode )
    {
        // Initialize output
        moris::Matrix< moris::IndexMat > tExpectedEdgeToNode = {
            { aElementToNode( iElem, 0 ), aElementToNode( iElem, 1 ) },
            { aElementToNode( iElem, 1 ), aElementToNode( iElem, 2 ) },
            { aElementToNode( iElem, 0 ), aElementToNode( iElem, 2 ) }
        };

        return tExpectedEdgeToNode;
    }

    //------------------------------------------------------------------------------

    inline bool
    verify_tri3_edge_topology(
            moris::Matrix< moris::IndexMat > const & aElementToNode,
            moris::Matrix< moris::IndexMat > const & aElementToEdge,
            moris::Matrix< moris::IndexMat > const & aEdgeToNode )
    {
        // Number of elements
        bool tValidEdgeTopo = true;

        moris::moris_index tNumElem        = aElementToNode.n_rows();
        moris::moris_index tNumEdgePerElem = aElementToEdge.n_cols();

        // Initilize expected edge to node connectivity
        moris::Matrix< moris::IndexMat > tExpectedEdgeToNode;
        moris::Matrix< moris::IndexMat > tActualEdgeToNode( 6, 2 );
        moris::Matrix< moris::IndexMat > tReorderEdgesMatrix( tNumElem, 6 );

        MORIS_ASSERT( tNumEdgePerElem == 3, "TRI3 NEEDS TO HAVE 3 EDGES (3 COLUMNS)" );

        for ( moris::moris_index i = 0; i < tNumElem; i++ )
        {
            tExpectedEdgeToNode = construct_expected_edge_to_node_tri3( i, aElementToNode );

            for ( moris::moris_index j = 0; j < tNumEdgePerElem; j++ )
            {
                moris::moris_index tEdgeIndex = aElementToEdge( i, j );

                // Make sure the same node does not appear on this edge twice.
                MORIS_ASSERT( aEdgeToNode( tEdgeIndex, 0 ) != aEdgeToNode( tEdgeIndex, 1 ),
                        "Same node cannot appear twice on the same edge." );

                // Loop over nodes in an edge and check for equivalence
                for ( moris::moris_index k = 0; k < 2; k++ )
                {
                    // Add to actual edge to node
                    tActualEdgeToNode( j, k ) = aEdgeToNode( tEdgeIndex, k );

                    // Check both nodes on edge ( assuming node order on an edge does not matter )
                    // Do something if node is invalid
                    if ( tExpectedEdgeToNode( j, k ) != aEdgeToNode( tEdgeIndex, 0 )
                            && tExpectedEdgeToNode( j, k ) != aEdgeToNode( tEdgeIndex, 1 ) )
                    {
                        std::cout << "iElem = " << i << " iEdge = " << j << std::endl;
                        std::cout << "Expected Nodes = " << tExpectedEdgeToNode( j, 0 ) << " " << tExpectedEdgeToNode( j, 1 ) << std::endl;
                        std::cout << "Actual Nodes   = " << aEdgeToNode( tEdgeIndex, 0 ) << " " << aEdgeToNode( tEdgeIndex, 1 ) << std::endl;
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
    verify_tri3_topology(
            moris::Matrix< moris::IndexMat > const & aElementToNode,
            moris::Matrix< moris::IndexMat > const & aElementToEdge,
            moris::Matrix< moris::IndexMat > const & aEdgeToNode )
    {
        MORIS_ASSERT( aElementToNode.n_cols() == 3, "UNDEFINED NUMBER OF NODES PROVIDED, NEEDS TO BE 3 FOR TRI3" );

        bool tValidTopo = false;

        // verify the edges
        bool tValidEdgeTopo = verify_tri3_edge_topology( aElementToNode, aElementToEdge, aEdgeToNode );

        if ( !tValidEdgeTopo )
        {
            std::cout << "Invalid edge topology detected" << std::endl;
            tValidTopo = false;
        }
        else
        {
            tValidTopo = true;
        }
        return tValidTopo;
    }

}    // namespace xtk

#endif /* SRC_MESH_FN_VERIFY_TRI_TOPOLOGY_HPP_ */
