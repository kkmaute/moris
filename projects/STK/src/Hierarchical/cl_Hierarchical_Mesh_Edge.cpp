/*
 * cl_Hierarchical_Mesh_Edge.cpp
 *
 *  Created on: Dec 12, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Edge.hpp"
using namespace moris;

Mat<uint>
Hierarchical_Mesh_Edge::give_edge_nodes(
        uint const & aEdgeId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Edge function is only available for 2D and 3D");
    // Compute the level of the edge
    uint tEdgeLevel = mBaseEdge.give_edge_level( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    // Compute the position of the edge
    Mat<uint> tEdgePosition = mBaseEdge.give_edge_position( aEdgeId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tEdgeNodes( 2 , 1 , UINT_MAX );
    //Potions for the nodes
    Mat<uint> tIJKPosition( aModelDim , 1 , 0 );
    if( aModelDim == 2 )
    {
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tEdgePosition( 0 );
            tIJKPosition( 1 ) = tEdgePosition( 1 );
            // Compute the edge nodes with the position of the edge and the functionality of the basis position, Corner nodes are only in linear possible
            tEdgeNodes( 0 ) = mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
            tIJKPosition( 1 ) = tEdgePosition( 1 );
            tEdgeNodes( 1 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tEdgePosition( 2 ) < UINT_MAX && tEdgePosition( 3 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tEdgePosition( 2 );
            tIJKPosition( 1 ) = tEdgePosition( 3 );
            tEdgeNodes( 0 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tEdgePosition( 2 );
            tIJKPosition( 1 ) = tEdgePosition( 3 ) + 1;
            tEdgeNodes( 1 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
        }
    }
    else if( aModelDim == 3 )
    {
        if( tEdgePosition( 0 ) < UINT_MAX && tEdgePosition( 1 ) < UINT_MAX && tEdgePosition( 2 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tEdgePosition( 0 );
            tIJKPosition( 1 ) = tEdgePosition( 1 );
            tIJKPosition( 2 ) = tEdgePosition( 2 );
            tEdgeNodes( 0 ) = mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tEdgePosition( 0 ) + 1;
            tIJKPosition( 1 ) = tEdgePosition( 1 );
            tIJKPosition( 2 ) = tEdgePosition( 2 );
            tEdgeNodes( 1 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tEdgePosition( 3 ) < UINT_MAX && tEdgePosition( 4 ) < UINT_MAX && tEdgePosition( 5 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tEdgePosition( 3 );
            tIJKPosition( 1 ) = tEdgePosition( 4 );
            tIJKPosition( 2 ) = tEdgePosition( 5 );
            tEdgeNodes( 0 ) = mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tEdgePosition( 3 );
            tIJKPosition( 1 ) = tEdgePosition( 4 ) + 1;
            tIJKPosition( 2 ) = tEdgePosition( 5 );
            tEdgeNodes( 1 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tEdgePosition( 6 ) < UINT_MAX && tEdgePosition( 7 ) < UINT_MAX && tEdgePosition( 8 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tEdgePosition( 6 );
            tIJKPosition( 1 ) = tEdgePosition( 7 );
            tIJKPosition( 2 ) = tEdgePosition( 8 );
            tEdgeNodes( 0 ) = mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tEdgePosition( 6 );
            tIJKPosition( 1 ) = tEdgePosition( 7 );
            tIJKPosition( 2 ) = tEdgePosition( 8 ) + 1;
            tEdgeNodes( 1 ) =  mBasis.give_basis_of_position( tEdgeLevel, aModelDim, 1, aNumberOfElementsPerDirection, tIJKPosition );
        }
    }
    return tEdgeNodes;
}
