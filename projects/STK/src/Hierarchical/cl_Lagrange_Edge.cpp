/*
 * cl_Lagrange_Edge.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#include "cl_Lagrange_Edge.hpp"
using namespace moris;

Mat<uint>
Lagrange_Edge::give_edge_nodes(
        uint const & aEdgeId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // Compute the position of the edge
    Mat<uint> tEdgePosition = mBaseEdge.give_edge_position( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the Elements, which are connected to the edge
    Mat<uint> tElementsOfEdge = mBaseEdge.give_elements_of_edge( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the nodes of the first element
    Mat<uint> tNodesOfElement = mLagrangeElement.give_basis_of_element( tElementsOfEdge( 0 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    Mat<uint> tEdgeNodes( aPolynomial + 1, 1, 0);
    if( aModelDim == 2 )
    {
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
        {
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tEdgeNodes( i ) = tNodesOfElement( tNodesOfElement.length() - aPolynomial - 1 + i );
            }
        }
        else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
        {
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tEdgeNodes( i ) = tNodesOfElement( aPolynomial + i * ( aPolynomial + 1 ) );
            }
        }
    }
    else if( aModelDim == 3 )
    {
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX && tEdgePosition( 2 ) < UINT_MAX )
        {
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tEdgeNodes( i ) = tNodesOfElement( tNodesOfElement.length() - aPolynomial - 1 + i );
            }
        }
        else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX )
        {
            //Compute nodes of a slice
            uint tNodesPerSlice = pow( aPolynomial + 1, 2 );
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tEdgeNodes( i ) = tNodesOfElement( ( 1 + i ) * tNodesPerSlice - 1 );
            }
        }
        else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX )
        {
            //Compute the first node of the z-direction
            uint tNodesToLastSlice = pow( aPolynomial + 1, 2 ) * aPolynomial + aPolynomial;
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tEdgeNodes( i ) = tNodesOfElement( tNodesToLastSlice + i * ( aPolynomial + 1 ) );
            }
        }
    }
    return tEdgeNodes;
}

