/*
 * cl_Lagrange_Face.cpp
 *
 *  Created on: Mar 1, 2018
 *      Author: gleim
 */

#include "cl_Lagrange_Face.hpp"
using namespace moris;
Mat<uint>
Lagrange_Face::give_face_basis(
        uint const & aFaceId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    MORIS_ASSERT( aModelDim > 1, "Face function is only available for 2D and 3D");
    // Compute the position of the Face
    Mat<uint> tFacePosition = mBaseFace.give_face_position( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the Elements, which are connected to the Face
    Mat<uint> tElementsOfFace = mBaseFace.give_elements_of_face( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the nodes of the first element
    Mat<uint> tNodesOfElement = mLagrangeElement.give_basis_of_element( tElementsOfFace( 0 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    Mat<uint> tFaceNodes;
    if( aModelDim == 2 )
    {
        tFaceNodes.set_size( aPolynomial + 1, 1, 0);
        //If face is in x-direction
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
        {
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tFaceNodes( i ) = tNodesOfElement( aPolynomial + i * ( aPolynomial + 1 ) );
            }
        }
        else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
        {
            //If face is in y-direction
            for ( uint i = 0; i <= aPolynomial; i++ )
            {
                tFaceNodes( i ) = tNodesOfElement( tNodesOfElement.length() - aPolynomial - 1 + i );
            }
        }
    }
    else if( aModelDim == 3 )
    {
        tFaceNodes.set_size( pow( aPolynomial + 1, 2 ), 1, 0);
        //Temporary variable for loop
        uint tVar = 0;
        //Compute nodes of a slice
        uint tNodesPerSlice = pow( aPolynomial + 1, 2 );
        //Compute nodes per edge
        uint tNodesPerEdge = aPolynomial + 1;
        //If face is in x-direction
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX && tFacePosition( 2 ) < UINT_MAX )
        {
            for ( uint j = 0; j <= aPolynomial; j++ )
            {
                for ( uint i = 0; i <= aPolynomial; i++ )
                {
                    tFaceNodes( tVar ) = tNodesOfElement( aPolynomial + tNodesPerSlice * j + tNodesPerEdge * i );
                    tVar++;
                }
            }
        }
        else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX )
        {
            //If face is in y-direction
            //Compute the number of nodes per edge until the last row of a slice
            uint tMultiTimesNodeEdges = ( aPolynomial + 1 ) * aPolynomial;
            for ( uint j = 0; j <= aPolynomial; j++ )
            {
                for ( uint i = 0; i <= aPolynomial; i++ )
                {
                    tFaceNodes( tVar ) = tNodesOfElement( tMultiTimesNodeEdges + i + j * tNodesPerSlice );
                    tVar++;
                }
            }
        }
        else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX )
        {
            //If face is in z-direction
            //Compute the first node of the last slize of nodes
            uint tNodesToLastSlice = pow( aPolynomial + 1, 2 ) * aPolynomial;
            for ( uint j = 0; j <= aPolynomial; j++ )
            {
                for ( uint i = 0; i <= aPolynomial; i++ )
                {
                    tFaceNodes( tVar ) = tNodesOfElement( tNodesToLastSlice + i + j * tNodesPerEdge );
                    tVar++;
                }
            }

        }
    }
    return tFaceNodes;
}

