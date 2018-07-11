/*
 * cl_Hierarchical_Mesh_Face.cpp
 *
 *  Created on: Dec 12, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Face.hpp"
using namespace moris;
Mat<uint>
Hierarchical_Mesh_Face::give_face_basis(
        uint const & aFaceId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    //Compute the level of the face id
    uint tFaceLevel = mBaseFace.give_face_level( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    //Compute the position of the face ID
    Mat<uint> tFacePosition = mBaseFace.give_face_position( aFaceId, aModelDim, aNumberOfElementsPerDirection );
    Mat<uint> tFaceBasis;
    Mat<uint> tIJKPosition( aModelDim, 1 , 0 );
    // Only corner nodes are needed, thats why aPolynomial is fixed to one!!!
    uint aPolynomial = 1;
    if( aModelDim == 2 )
    {
        tFaceBasis.set_size( 2, 1, UINT_MAX );
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tFacePosition( 0 );
            tIJKPosition( 1 ) = tFacePosition( 1 );
            tFaceBasis( 0 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 0 );
tIJKPosition( 1 ) = tFacePosition( 1 ) + 1;
            tFaceBasis( 1 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tFacePosition( 2 ) < UINT_MAX && tFacePosition( 3 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tFacePosition( 2 );
            tIJKPosition( 1 ) = tFacePosition( 3 );
            tFaceBasis( 0 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 2 ) + 1;
            tIJKPosition( 1 ) = tFacePosition( 3 );
            tFaceBasis( 1 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
        }
    }
    else if( aModelDim == 3 )
    {
        tFaceBasis.set_size( 4, 1, UINT_MAX );
        if( tFacePosition( 0 ) < UINT_MAX && tFacePosition( 1 ) < UINT_MAX && tFacePosition( 2 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tFacePosition( 0 );
            tIJKPosition( 1 ) = tFacePosition( 1 );
            tIJKPosition( 2 ) = tFacePosition( 2 );
            tFaceBasis( 0 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 0 );
            tIJKPosition( 1 ) = tFacePosition( 1 ) + 1;
            tIJKPosition( 2 ) = tFacePosition( 2 );
            tFaceBasis( 1 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 0 );
            tIJKPosition( 1 ) = tFacePosition( 1 );
            tIJKPosition( 2 ) = tFacePosition( 2 ) + 1;
            tFaceBasis( 2 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 0 );
            tIJKPosition( 1 ) = tFacePosition( 1 ) + 1;
            tIJKPosition( 2 ) = tFacePosition( 2 ) + 1;
            tFaceBasis( 3 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tFacePosition( 3 ) < UINT_MAX && tFacePosition( 4 ) < UINT_MAX && tFacePosition( 5 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tFacePosition( 3 );
            tIJKPosition( 1 ) = tFacePosition( 4 );
            tIJKPosition( 2 ) = tFacePosition( 5 );
            tFaceBasis( 0 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 3 ) + 1;
            tIJKPosition( 1 ) = tFacePosition( 4 );
            tIJKPosition( 2 ) = tFacePosition( 5 );
            tFaceBasis( 1 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 3 );
            tIJKPosition( 1 ) = tFacePosition( 4 );
            tIJKPosition( 2 ) = tFacePosition( 5 ) + 1;
            tFaceBasis( 2 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 3 ) + 1;
            tIJKPosition( 1 ) = tFacePosition( 4 );
            tIJKPosition( 2 ) = tFacePosition( 5 ) + 1;
            tFaceBasis( 3 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
        }
        else if( tFacePosition( 6 ) < UINT_MAX && tFacePosition( 7 ) < UINT_MAX && tFacePosition( 8 ) < UINT_MAX )
        {
            tIJKPosition( 0 ) = tFacePosition( 6 );
            tIJKPosition( 1 ) = tFacePosition( 7 );
            tIJKPosition( 2 ) = tFacePosition( 8 );
            tFaceBasis( 0 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 6 ) + 1;
            tIJKPosition( 1 ) = tFacePosition( 7 );
            tIJKPosition( 2 ) = tFacePosition( 8 );
            tFaceBasis( 1 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 6 );
            tIJKPosition( 1 ) = tFacePosition( 7 ) + 1;
            tIJKPosition( 2 ) = tFacePosition( 8 );
            tFaceBasis( 2 ) = mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            tIJKPosition( 0 ) = tFacePosition( 6 ) + 1;
            tIJKPosition( 1 ) = tFacePosition( 7 ) + 1;
            tIJKPosition( 2 ) = tFacePosition( 8 );
            tFaceBasis( 3 ) =  mBasis.give_basis_of_position( tFaceLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
        }
    }
    return tFaceBasis;
}
