/*
 * cl_Lagrange_Basis.cpp
 *
 *  Created on: Feb 21, 2018
 *      Author: gleim
 */

#include "cl_Lagrange_Basis.hpp"
using namespace moris;

/*
uint
Lagrange_Basis::give_number_of_basis(
        uint const & aLevel,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection)
{
    // Variable to calculate the number of basis functions from level 0 until level "aLevel"
    uint tBasisNumber = 0;
    if( aModelDim == 1 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tBasisNumber += pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 ;
        }
    }
    else if( aModelDim == 2 )
    {
        for( uint i = 0; i < aLevel + 1; i++ )
        {
            tBasisNumber +=  ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                            * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
        }
    }
    else if( aModelDim == 3 )
    {
        for( uint i = 0; i < aLevel + 1 ; i++ )
        {
            tBasisNumber += ( pow( 2, i ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                            * ( pow( 2, i ) * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                                            * ( pow( 2, i ) * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
        }
    }
    return tBasisNumber;
} */

//--------------------------------------------------------------------------------

Mat<uint>
Lagrange_Basis::give_neighbor_of_basis(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        uint const & aBuffer,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    //Position of the basis function
    Mat<uint> tBasisPosition = give_position_of_basis( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    // Determine the level of the basis function
    uint tBasisLevel = give_basis_level( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    // Element plus neighbors (9 for aModelDim = 2, 27 for aModelDim = 3 for a buffer of 1, otherwise more)
    Mat<uint> tBasisNeighbor( pow( 1 + 2 * aBuffer, aModelDim ) , 1 , UINT_MAX );
    // Position of a neighbour element;
    Mat<uint> tNeighborPosition( aModelDim, 1 , 0 );
    uint tPowVar = pow( 2, tBasisLevel );
    //Compute the Lower and Upper bound of possible basis functions
    uint tBasisLowerBoundPerDirection = aPolynomial * tPowVar;
    Mat<uint> tBasisUpperBoundPerDirection( aModelDim, 1, 0 );
    if( aModelDim == 1 )
    {
        tBasisUpperBoundPerDirection(0) = aNumberOfElementsPerDirection(0) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        if( aBuffer == 1 )
        {
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tBasisNeighbor( 0 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tBasisNeighbor( 1 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tBasisNeighbor( 2 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else
        {
            // If a higher buffer layer is needed, a loop generates the neighbours
            uint tVar = 0;
            Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
            for( int i = -( (int)aBuffer ); i <= ( (int)aBuffer ); i++ )
            {
                if ( tBasisPosition( 0 ) + i >= tBasisLowerBoundPerDirection
                        && tBasisPosition( 0 ) + i < tBasisUpperBoundPerDirection( 0 ) )
                {
                    tNeighborPosition( 0 ) = tBasisPosition( 0 ) + i;
                    tBasisNeighbor(tVar) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
                }
                tVar++;
            }
        }
    }
    else if( aModelDim == 2 )
    {
        tBasisUpperBoundPerDirection(0) = aNumberOfElementsPerDirection(0) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        tBasisUpperBoundPerDirection(1) = aNumberOfElementsPerDirection(1) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        if( aBuffer == 1 )
        {
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tBasisNeighbor( 0 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tBasisNeighbor( 1 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tBasisNeighbor( 2 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tBasisNeighbor( 3 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tBasisNeighbor( 4 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tBasisNeighbor( 5 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tBasisNeighbor( 6 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tBasisNeighbor( 7 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tBasisNeighbor( 8 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else
        {
            // If a higher buffer layer is needed, a loop generates the neighbours
            uint tVar = 0;
            Mat<uint> aIJKPosition(aModelDim,1,0);
            for( int j = -( (int)aBuffer ); j <= ( (int)aBuffer ); j++ )
            {
                for(int i = -( (int)aBuffer ); i <= ( (int)aBuffer ); i++ )
                {
                    if ( tBasisPosition( 0 ) + i >= tBasisLowerBoundPerDirection
                            && tBasisPosition( 0 ) + i <  tBasisUpperBoundPerDirection( 0 )
                            && tBasisPosition( 1 ) + j >= tBasisLowerBoundPerDirection
                            && tBasisPosition( 1 ) + j <  tBasisUpperBoundPerDirection( 1 ) )
                    {
                        tNeighborPosition( 0 ) = tBasisPosition( 0 ) + i;
                        tNeighborPosition( 1 ) = tBasisPosition( 1 ) + j;
                        tBasisNeighbor(tVar) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
                    }
                    tVar++;
                }
            }
        }
    }
    else if( aModelDim == 3 )
    {
        tBasisUpperBoundPerDirection(0) = aNumberOfElementsPerDirection(0) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        tBasisUpperBoundPerDirection(1) = aNumberOfElementsPerDirection(1) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        tBasisUpperBoundPerDirection(2) = aNumberOfElementsPerDirection(2) * tPowVar + aPolynomial - aPolynomial * tPowVar;
        //If buffer == 1, find the direct neighbor elements
        if( aBuffer == 1 )
        {
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 0 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 1 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 2 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 3 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 4 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 5 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 6 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 7 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) - 1;
            tBasisNeighbor( 8 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 9 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 10 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 11 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 12 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 13 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 14 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 15 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 16 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 );
            tBasisNeighbor( 17 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 18 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 19 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 20 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 21 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 22 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 );
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 23 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 24 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 );
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 25 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + 1;
            tBasisNeighbor( 26 ) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else
        {
            // If a higher buffer layer is needed, a loop generates the neighbors
            uint tVar = 0;
            Mat<uint> aIJKPosition( aModelDim, 1, 0 );
            for( int k = -( (int)aBuffer ); k <= ( (int)aBuffer ); k++ )
            {
                for( int j = -( (int)aBuffer ); j <= ( (int)aBuffer ); j++ )
                {
                    for( int i = -( (int)aBuffer ); i <= ( (int)aBuffer ); i++ )
                    {
                        if ( tBasisPosition( 0 )+i >= tBasisLowerBoundPerDirection
                                && tBasisPosition( 0 ) + i < tBasisUpperBoundPerDirection( 0 )
                                && tBasisPosition( 1 )+j >= tBasisLowerBoundPerDirection
                                && tBasisPosition( 1 ) + j < tBasisUpperBoundPerDirection( 1 )
                                && tBasisPosition( 2 )+k >= tBasisLowerBoundPerDirection
                                && tBasisPosition( 2 ) + k < tBasisUpperBoundPerDirection( 2 ) )
                        {
                            tNeighborPosition( 0 ) = tBasisPosition( 0 ) + i;
                            tNeighborPosition( 1 ) = tBasisPosition( 1 ) + j;
                            tNeighborPosition( 2 ) = tBasisPosition( 2 ) + k;
                            tBasisNeighbor(tVar) = give_basis_of_position( tBasisLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tNeighborPosition );
                        }
                        tVar++;
                    }
                }
            }
        }
    }
    return tBasisNeighbor;
}

//--------------------------------------------------------------------------------

/* uint
Lagrange_Basis::give_basis_level(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tBasisLevel = 1;
    // temporary variable for the while loop
    uint tLevel = 1;
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = pow( 2, tBasisLevel - 1 );
    if( aModelDim == 1 )
    {
        // temporary variable for the number of basis per level
        uint tNumberOfBasis = tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
        while(tLevel>0)
        {
            // Check if it fits in the level
            tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
            if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                tBasisLevel++;
                // Count the number of basis functions on the current level
                tPowLevel = pow( 2, tBasisLevel - 1 );
                tNumberOfBasis += tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
            }
        }
        tBasisLevel--;
    }
    else if( aModelDim == 2 )
    {
        // temporary variable for the number of basis per level
        uint tNumberOfBasis = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                            * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
        while( tLevel > 0 )
        {
            // Check if it fits in the level
            tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
            if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                tBasisLevel++;
                // Count the number of basis functions on the current level
                tPowLevel = pow( 2, tBasisLevel - 1 );
                tNumberOfBasis += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
            }
        }
        tBasisLevel--;
    }
    else if( aModelDim == 3 )
    {
        // temporary variable for the number of elements per level
        uint tNumberOfBasis = ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                            * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                            * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
        while( tLevel > 0 )
        {
            // Check if it fits in the level
            tLevel = floor( (real)( aBasisId + 1 ) / tNumberOfBasis );
            if( (real)( aBasisId + 1 ) / tNumberOfBasis == 1 || tLevel < 1 )
            {
                tLevel = 0;
            }
            else
            {
                tBasisLevel++;
                // Count the number of basis functions on the current level
                tPowLevel = pow( 2, tBasisLevel - 1 );
                tNumberOfBasis += ( tPowLevel * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                                                                                                * ( tPowLevel * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 )
                                                                                                * ( tPowLevel * aNumberOfElementsPerDirection( 2 ) * aPolynomial + 1 );
            }
        }
        tBasisLevel--;
    }
    return tBasisLevel;
} */

//--------------------------------------------------------------------------------

/*
Mat<uint>
Lagrange_Basis::give_position_of_basis(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Basis positions in the direction x,y,z
    Mat<uint> tBasisPosition( aModelDim , 1 , 0 );
    // Compute the level of the basis function
    uint tBasisLevel=give_basis_level( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    // temporary variable for the number of elements per level
    uint tNumberOfBasis = 0;
    if( aModelDim == 2 )
    {
        if( tBasisLevel != 0 )
        {
            //Count basis functions until tBasisLevel-1
            tNumberOfBasis = give_number_of_basis(tBasisLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        }

        tBasisPosition( 1 ) = ceil( (real)(aBasisId + 1 - tNumberOfBasis )
                / (real)( ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 ) ) );
        // Calculate the x-direction with the calculated y-position
        tBasisPosition( 0 ) = aBasisId +1 - tNumberOfBasis
                + ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 )
                - tBasisPosition( 1 ) * ( pow( 2, tBasisLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 );
        // -1 to get an indexed basis position, starting with zero
        tBasisPosition( 1 ) -= 1;
        tBasisPosition( 0 ) -= 1;
    }
    else if( aModelDim == 3 )
    {
        if( tBasisLevel !=0 )
        {
            //Count basis functions until tBasisLevel-1
            tNumberOfBasis = give_number_of_basis( tBasisLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        }
        //Temporary variables with precalculation for simpler coding style
        uint tPow2 = pow( 2, tBasisLevel );

        uint tConst0 = tPow2 * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
        uint tConst1 = tPow2 * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1;
        // Calculate the position in z-direction
        tBasisPosition( 2 ) = ceil( (real)( aBasisId + 1 - tNumberOfBasis ) / ( (real)( tConst0 * tConst1 ) ) );
        // Calculate the y-direction with the calculated z-position
        tBasisPosition( 1 ) = ceil( (real)( aBasisId + 1 - tNumberOfBasis - ( tBasisPosition( 2 ) - 1 ) * tConst0 * tConst1 ) / (real)tConst0 );
        // Calculate the x-direction with the calculated y and z-position
        tBasisPosition( 0 ) = aBasisId + 1 - tNumberOfBasis + tConst0 * ( 1 - tBasisPosition( 1 ) + tConst1 * ( 1 - tBasisPosition( 2 ) ) );
        // -1 to get an indexed basis position, starting with zero
        tBasisPosition( 2 ) -= 1;
        tBasisPosition( 1 ) -= 1;
        tBasisPosition( 0 ) -= 1;
    }
    return tBasisPosition;
} */

//--------------------------------------------------------------------------------

/*
uint
Lagrange_Basis::give_basis_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    uint tBasisId = 0;
    // Position of the Element in y-direction
    uint tY_Position = 0;
    if( aModelDim == 1 )
    {
        if( aLevel > 0 )
        {
            //Give number of basis functions for the last level to have the initial number at position 0,0
            tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        }
        //Basis function with the number of basis of lower level
        tBasisId += aIJKPosition( 0 );
    }
    else if( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 )*( aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1 );
        if( aLevel>0 )
        {
            //Give number of basis functions for the last level to have the initial number at position 0,0
            tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
            // Position of the Element in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( 1 + aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) );
        }
        //Basis function with the calculated y-position
        tBasisId += aIJKPosition( 0 ) + tY_Position;
    }
    else if( aModelDim == 3 )
    {
        tY_Position = aIJKPosition( 1 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) + 1 )
                                                                                           + aIJKPosition( 2 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) + 1 )
                                                                                           * ( aPolynomial * aNumberOfElementsPerDirection( 1 ) + 1 );
        if( aLevel > 0 )
        {
            //Give number of basis functions for the last level to have the initial number at position 0,0,0
            tBasisId = give_number_of_basis( aLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
            // Position of the Element in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) +1  )
                                                                                                            + aIJKPosition( 2 ) * ( aPolynomial * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel ) + 1 )
                                                                                                            * ( aPolynomial * aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel ) + 1 );
        }
        tBasisId += aIJKPosition( 0 ) + tY_Position; //Basis function with the calculated y-position
    }
    return tBasisId;
} */

//--------------------------------------------------------------------------------

uint
Lagrange_Basis::give_basis_of_parent(
        uint const & aBasis,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    Mat<uint> tIJKPosition = give_position_of_basis( aBasis, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    uint tBasis = UINT_MAX;
    uint tLevel = give_basis_level( aBasis, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    if( tLevel > 0 )
    {
        // Parent basis function is on a lower level!
        tLevel -= 1;
        if( aModelDim == 2 )
        {
            // Check if a Position has a even number, otherwise there is no parent existent
            if( tIJKPosition( 0 ) % 2 == 0 && tIJKPosition( 1 ) % 2 == 0 )
            {
                // Divide by two, to get the position on the coarser level
                tIJKPosition( 0 ) /= 2;
                tIJKPosition( 1 ) /= 2;
                // Compute the basis id with the position
                tBasis = give_basis_of_position( tLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            }
        }
        else if( aModelDim == 3 )
        {
            // Check if a Position has a even number, otherwise there is no parent existent
            if( tIJKPosition( 0 ) % 2 == 0 && tIJKPosition( 1 ) % 2 == 0 && tIJKPosition( 2 ) % 2 == 0 )
            {
                // Divide by two, to get the position on the coarser level
                tIJKPosition( 0 ) /= 2;
                tIJKPosition( 1 ) /= 2;
                tIJKPosition( 2 ) /= 2;
                // Compute the basis id with the position
                tBasis = give_basis_of_position( tLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection, tIJKPosition );
            }
        }
    }
    return tBasis;
}

//--------------------------------------------------------------------------------

Mat<uint>
Lagrange_Basis::give_element_of_basis(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Compute level of basis function
    uint tBasisLevel = give_basis_level( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    // Compute Position of basis function
    Mat<uint> tBasisPosition = give_position_of_basis( aBasisId, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    // Elements, which have support with the basis function aBasis
    Mat<uint> tElements( pow( 2, aModelDim ) , 1, UINT_MAX );
    // Position of the element with respect to the basis function
    Mat<uint> tElementPosition( aModelDim , 1 , 0 );
    //Temporary variable
    uint tPowVar = pow( 2, tBasisLevel );
    //Temporary variable for if clause conditions, it is the maximum number of elements in a specific direction
    uint tDomainRangeCheck0 = tPowVar * ( aNumberOfElementsPerDirection( 0 ) ) * aPolynomial + 1;
    if( aModelDim == 1 )
    {
        // Elements are hardcoded for aPolynomial = 1 & 2. aPolynomial > 2 is generated by a loop
        // Temporary variable for the for loops
        uint tVar=0;
        for( uint i = 0; i < 2; i++ )
        {
            if( (int)( tBasisPosition( 0 ) + i ) > 0
                    && (int)( tBasisPosition( 0 ) + i ) < (int)tDomainRangeCheck0)
            {
                tElementPosition( 0 ) = floor( (real)( tBasisPosition(0) + i )
                        / (real)tDomainRangeCheck0 * (real)( aNumberOfElementsPerDirection( 0 ) * tPowVar ) );
                tElements( tVar ) = mBaseElement.give_element_of_position( tBasisLevel, aModelDim, aNumberOfElementsPerDirection, tElementPosition );
            }
            tVar++;
        }
    }
    else if( aModelDim == 2 )
    {
        //Temporary variable for if clause conditions, it is the maximum number of elements in a specific direction
        uint tDomainRangeCheck1 = tPowVar * ( aNumberOfElementsPerDirection( 1 ) ) * aPolynomial + 1;
        // Temporary variable for the for loops
        uint tVar=0;
        for( uint j = 0; j < 2; j++ )
        {
            for( uint i = 0; i < 2; i++ )
            {
                if( (int)( tBasisPosition( 0 ) + i ) > 0 && (int)( tBasisPosition( 1 ) + j ) > 0
                        && (int)( tBasisPosition( 0 ) + i ) < (int)tDomainRangeCheck0
                        && (int)( tBasisPosition( 1 ) + j ) < (int)tDomainRangeCheck1 )
                {
                    tElementPosition( 0 ) = floor( (real)( tBasisPosition(0) + i )
                            / (real)tDomainRangeCheck0 * (real)( aNumberOfElementsPerDirection( 0 ) * tPowVar ) );
                    tElementPosition( 1 ) = floor( (real)( tBasisPosition(1) + j )
                            / (real)tDomainRangeCheck1 * (real)( aNumberOfElementsPerDirection( 1 ) * tPowVar ) );
                    tElements( tVar ) = mBaseElement.give_element_of_position( tBasisLevel, aModelDim, aNumberOfElementsPerDirection, tElementPosition );
                }
                tVar++;
            }
        }

    }
    else if( aModelDim == 3 )
    {
        //Temporary variable for if clause conditions, it is the maximum number of elements in a specific direction
        uint tDomainRangeCheck1 = tPowVar * ( aNumberOfElementsPerDirection( 1 ) ) * aPolynomial + 1;
        uint tDomainRangeCheck2 = tPowVar * ( aNumberOfElementsPerDirection( 2 ) ) * aPolynomial + 1;
        // Elements are hardcoded for aPolynomial = 1. aPolynomial > 1 is generated by a loop
        uint tVar = 0; // Temporary variable for the for loops
        for ( uint k = 0; k < 2; k++)
        {
            for( uint j = 0; j < 2; j++)
            {
                for ( uint i = 0; i < 2; i++)
                {
                    if((int)(tBasisPosition( 0 ) + i) > 0 && (int)(tBasisPosition( 1 ) + j) > 0
                            && (int)(tBasisPosition( 2 ) + k) > 0
                            && (int)(tBasisPosition( 0 ) + i) < (int)tDomainRangeCheck0
                            && (int)(tBasisPosition( 1 ) + j) < (int)tDomainRangeCheck1
                            && (int)(tBasisPosition( 2 ) + k) < (int)tDomainRangeCheck2 )
                    {
                        tElementPosition( 0 ) = floor( (real)( tBasisPosition(0) + i )
                                / (real)tDomainRangeCheck0 * (real)( aNumberOfElementsPerDirection( 0 ) * tPowVar ) );
                        tElementPosition( 1 ) = floor( (real)( tBasisPosition(1) + j )
                                / (real)tDomainRangeCheck1 * (real)( aNumberOfElementsPerDirection( 1 ) * tPowVar ) );
                        tElementPosition( 2 ) = floor( (real)( tBasisPosition(2) + j )
                                / (real)tDomainRangeCheck2 * (real)( aNumberOfElementsPerDirection( 2 ) * tPowVar ) );
                        (tElements(tVar)) = mBaseElement.give_element_of_position( tBasisLevel, aModelDim, aNumberOfElementsPerDirection, tElementPosition );
                    }
                    tVar++;
                }
            }
        }
    }
    return tElements;
}

//--------------------------------------------------------------------------------

Mat<real>
Lagrange_Basis::give_coordinate_from_basis(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<real> const & aDomainRange,
        Mat<real> const & aDomainRangeOffset) const
{
    Mat<real> tCoordinates(aModelDim,1,0);
    // Compute the level of the basis function
    uint tBasisLevel = give_basis_level(aBasisId,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
    // Determine the position of the basis function
    Mat<uint> tBasisPosition = give_position_of_basis(aBasisId,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
    //Calculate the Coordinates of the basis function with respect to the dimensions and the offset
    //Temporary variables for precalculation
    real tPow2Level = pow(2,tBasisLevel);
    for( uint k = 0; k < aModelDim; k++ )
    {
        tCoordinates(k) = aDomainRangeOffset(k) + tBasisPosition(k) * aDomainRange(k) / ( tPow2Level * aNumberOfElementsPerDirection(k) * aPolynomial );
    }
    return tCoordinates;
}

//--------------------------------------------------------------------------------

uint
Lagrange_Basis::give_basis_proc_owner(
        uint const & aBasisId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcNeighbor,
        Mat<uint> const & aDecomp) const
{
    uint tProcRank = par_rank();
    uint tBasisLevel;
    Mat<uint> tBasisPosition( aModelDim , 1 , 0 );
    uint tBasisRank = UINT_MAX;
    if ( par_size() > 1)
    {
        // Compute the level of the basis function
        tBasisLevel = give_basis_level(aBasisId,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        // Compute the position of the basis function
        tBasisPosition = give_position_of_basis(aBasisId,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        Mat<uint> tProcRange(aDecomp.length(),1,0);
        // Compute the range on the respective level
        tProcRange( 0 ) = ( aDecomp( 0 ) - 1 ) * pow( 2, tBasisLevel ) + 1;
        tProcRange( 1 ) = aDecomp( 1 ) * pow( 2, tBasisLevel );
        tProcRange( 2 ) = ( aDecomp( 2 ) - 1 ) * pow( 2, tBasisLevel ) + 1;
        tProcRange( 3 ) = aDecomp( 3 ) * pow( 2, tBasisLevel );
        //See the numbering and the position of the processors in Hierarchical_Mesh_Main::proc_coordinates_neighbors
        if( aModelDim == 2 )
        {
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) >= tProcRange( 2 ) && tBasisPosition( 1 ) <= tProcRange( 3 ) ) // Inner basis
            {
                tBasisRank = tProcRank;
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) && tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) && aProcNeighbor( 3 ) <UINT_MAX ) // Left edge
            {
                tBasisRank = aProcNeighbor( 3 );
            }
                        else if ( tBasisPosition( 0 ) > tProcRange( 1 ) && tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) && aProcNeighbor( 1 ) <UINT_MAX ) // Right edge
            {
                tBasisRank = aProcNeighbor( 1 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) < tProcRange( 2 ) && aProcNeighbor( 0 ) < UINT_MAX ) // Bottom
            {
                tBasisRank = aProcNeighbor( 0 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) > tProcRange( 3 ) && aProcNeighbor( 2 ) < UINT_MAX ) // Top
            {
                tBasisRank = aProcNeighbor( 2 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) && tBasisPosition( 1 ) < tProcRange( 2 )
                    && aProcNeighbor( 8 ) <UINT_MAX ) // Bottom left
            {
                tBasisRank = aProcNeighbor( 8 );
            }
                        else if ( tBasisPosition( 0 ) > tProcRange( 1 ) && tBasisPosition( 1 ) < tProcRange( 2 )
                    && aProcNeighbor( 15 ) <UINT_MAX ) // Bottom right
            {
                tBasisRank = aProcNeighbor( 15 );
            }
                        else if ( tBasisPosition( 0 ) < tProcRange( 0 ) && tBasisPosition( 1 ) > tProcRange( 3 )
                    && aProcNeighbor( 16 ) <UINT_MAX ) // Top left
            {
                tBasisRank = aProcNeighbor( 16 );
            }
                                    else if ( tBasisPosition( 0 ) > tProcRange( 1 ) && tBasisPosition( 1 ) > tProcRange( 3 )
                    && aProcNeighbor( 17 ) <UINT_MAX ) // Top right
            {
                tBasisRank = aProcNeighbor( 17 );
            }
          }
        else if( aModelDim == 3)
        {
            // Compute the range on the respective level
            tProcRange( 4 ) = ( aDecomp( 4 ) - 1 ) * pow( 2, tBasisLevel ) + 1;
            tProcRange( 5 ) = aDecomp( 5 ) * pow( 2, tBasisLevel );
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) >= tProcRange( 4 ) && tBasisPosition( 2 ) <= tProcRange( 5 )) // Inner basis
            {
                tBasisRank = tProcRank;
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                    &&  tBasisPosition( 2 ) >= tProcRange( 4 ) && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 0 ) <UINT_MAX ) // Bottom
            {
                tBasisRank = aProcNeighbor( 0 );
            }
                        else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 )
                    &&  tBasisPosition( 2 ) >= tProcRange( 4 ) && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 2 ) <UINT_MAX ) // Top
            {
                tBasisRank = aProcNeighbor( 2 );
            }
                       else if ( tBasisPosition( 0 ) < tProcRange( 0 )  &&  tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) >= tProcRange( 4 ) && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 3 ) <UINT_MAX ) // Left
            {
                tBasisRank = aProcNeighbor( 3 );
            }
                                   else if ( tBasisPosition( 0 ) > tProcRange( 1 )  &&  tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) >= tProcRange( 4 ) && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 1 ) <UINT_MAX ) // Right
            {
                tBasisRank = aProcNeighbor( 1 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) < tProcRange( 4 ) && aProcNeighbor( 4 ) <UINT_MAX ) // Back
            {
                tBasisRank = aProcNeighbor( 4 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 )
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) > tProcRange( 5 ) && aProcNeighbor( 5 ) <UINT_MAX ) // Front
            {
                tBasisRank = aProcNeighbor( 5 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 6 ) <UINT_MAX ) // Back bottom left
            {
                tBasisRank = aProcNeighbor( 6 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 9 ) <UINT_MAX ) // Back bottom
            {
                tBasisRank = aProcNeighbor( 9 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 10 ) <UINT_MAX ) // Back bottom right
            {
                tBasisRank = aProcNeighbor( 10 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 7 ) <UINT_MAX ) // Back left
            {
                tBasisRank = aProcNeighbor( 7 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 11 ) <UINT_MAX ) // Back right
            {
                tBasisRank = aProcNeighbor( 11 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 12 ) <UINT_MAX ) // Back top left
            {
                tBasisRank = aProcNeighbor( 12 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 13 ) <UINT_MAX ) // Back top
            {
                tBasisRank = aProcNeighbor( 13 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) < tProcRange( 4 )  && aProcNeighbor( 14 ) <UINT_MAX ) // Back top right
            {
                tBasisRank = aProcNeighbor( 14 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) < tProcRange( 2 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 )  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 8 ) <UINT_MAX ) // Bottom left
            {
                tBasisRank = aProcNeighbor( 8 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 )  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 15 ) <UINT_MAX ) // Bottom right
            {
                tBasisRank = aProcNeighbor( 15 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 )  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 16 ) <UINT_MAX ) // Top left
            {
                tBasisRank = aProcNeighbor( 16 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 )  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 17 ) <UINT_MAX ) // Top right
            {
                tBasisRank = aProcNeighbor( 17 );
            }
                        else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 18 ) <UINT_MAX ) // Front bottom left
            {
                tBasisRank = aProcNeighbor( 18 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 19 ) <UINT_MAX ) // Front bottom
            {
                tBasisRank = aProcNeighbor( 19 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) < tProcRange( 2 )
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 20 ) <UINT_MAX ) // Front bottom right
            {
                tBasisRank = aProcNeighbor( 20 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 21 ) <UINT_MAX ) // Front left
            {
                tBasisRank = aProcNeighbor( 21 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 22 ) <UINT_MAX ) // Front right
            {
                tBasisRank = aProcNeighbor( 22 );
            }
            else if ( tBasisPosition( 0 ) < tProcRange( 0 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 23 ) <UINT_MAX ) // Front top left
            {
                tBasisRank = aProcNeighbor( 23 );
            }
            else if ( tBasisPosition( 0 ) >= tProcRange( 0 ) && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 24 ) <UINT_MAX ) // Front top
            {
                tBasisRank = aProcNeighbor( 24 );
            }
            else if ( tBasisPosition( 0 ) > tProcRange( 1 ) &&  tBasisPosition( 1 ) > tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) > tProcRange( 5 )  && aProcNeighbor( 25 ) <UINT_MAX ) // Front top right
            {
                tBasisRank = aProcNeighbor( 25 );
            }
         }
    }
    return tBasisRank;
}

//--------------------------------------------------------------------------------

Mat<uint>
Lagrange_Basis::give_basis_proc_owner(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aNodalLocaltoGlobal,
        Mat<uint> const & aProcNeighbour,
        Mat<uint> const & aDecomp,
        uint const & aNumberBasis,
        BoostBitset const & aBasisActive) const
{
    uint tProcRank = par_rank();
    Mat<uint> tNodeProcs( aNodalLocaltoGlobal.length() , 1 , tProcRank );
    if ( par_size() > 1 )
    {
        if ( tProcRank > 0 )
        {
            for( uint i = 0; i < aNodalLocaltoGlobal.length(); i++ )
            {
                if( aNodalLocaltoGlobal(i) < aBasisActive.size() && aBasisActive.test(aNodalLocaltoGlobal(i)) == 1 )
                {
                    tNodeProcs(i) = give_basis_proc_owner(aNodalLocaltoGlobal(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection,aProcNeighbour,aDecomp);
                }
            }
        }
    }
    return tNodeProcs;
}

//--------------------------------------------------------------------------------

Mat<uint>
Lagrange_Basis::give_basis_share(
        uint const & aBasis,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aProcNeighbor,
        Mat<uint> const& aDecomp) const
{
    uint tProcRank = par_rank();
    uint tBasisLevel;
    Mat<uint> tBasisPosition( aModelDim , 1 , 0 );
    uint tVar = 0;
    // Arbitrary number for possible sharing elements ( It will find with thie algorithm muliple times the same neighbors and unique will kick them out)
    Mat<uint> tProcShare( 30 , 1 , UINT_MAX );
    if ( par_size() > 1 )
    {
        // Compute the level of the basis function
        tBasisLevel = give_basis_level(aBasis,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        // Compute the position of the basis function
        tBasisPosition = give_position_of_basis(aBasis,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        // Compute the range on the respective level
        Mat<uint> tProcRange( aDecomp.length() , 1 , 0 );
        tProcRange( 0 ) = ( aDecomp( 0 ) - 1 ) * pow( 2, tBasisLevel )  + 1;
        tProcRange( 1 ) = aDecomp( 1 ) * pow( 2, tBasisLevel );
        tProcRange( 2 ) = ( aDecomp( 2 ) - 1 ) * pow( 2, tBasisLevel )  + 1;
        tProcRange( 3 ) = aDecomp( 3 ) * pow( 2, tBasisLevel );
        if( aModelDim == 2)
        {
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) >= tProcRange( 2 ) - 1 && tBasisPosition( 1 ) <= tProcRange( 3 ) ) // Inner basis
            {
                tProcShare(tVar) = tProcRank;
                 tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1  && tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) && aProcNeighbor( 3 ) <UINT_MAX ) //Left
            {
                tProcShare(tVar) = aProcNeighbor( 3 );
                 tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 )  && tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) && aProcNeighbor( 1 ) <UINT_MAX ) // Right
            {
                tProcShare(tVar) = aProcNeighbor( 1 );
                 tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) == tProcRange( 2 ) - 1 && aProcNeighbor( 0 ) <UINT_MAX ) // Bottom
            {
                tProcShare(tVar) = aProcNeighbor( 0 ); 
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 )
                    && tBasisPosition( 1 ) == tProcRange( 3 ) && aProcNeighbor( 2 ) <UINT_MAX ) // Top
            {
                tProcShare(tVar) = aProcNeighbor( 2 );
                 tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 && tBasisPosition( 1 ) == tProcRange( 2 ) - 1 && aProcNeighbor( 8 ) <UINT_MAX ) // Bottom left
            {
                tProcShare(tVar) = aProcNeighbor( 8 );
                 tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) && tBasisPosition( 1 ) == tProcRange( 2 ) - 1 && aProcNeighbor( 15 ) <UINT_MAX ) // Bottom right
            {
                tProcShare(tVar) = aProcNeighbor( 15 );
                 tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 && tBasisPosition( 1 ) == tProcRange( 3 ) && aProcNeighbor( 16 ) <UINT_MAX ) // Top left
            {
                tProcShare(tVar) = aProcNeighbor( 16 ); 
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) && tBasisPosition( 1 ) == tProcRange( 3 ) && aProcNeighbor( 17 ) <UINT_MAX ) // Top right
            {
                tProcShare(tVar) = aProcNeighbor( 17 ); 
                tVar++;
            }
                    }
        else if( aModelDim == 3 )
        {
            tProcRange( 4 ) = ( aDecomp( 4 ) - 1 ) * pow( 2, tBasisLevel ) + 1;
            tProcRange( 5 ) = aDecomp( 5 ) * pow( 2, tBasisLevel );
      if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1 && tBasisPosition( 2 ) <= tProcRange( 5 )) // Inner basis
            {
                tProcShare(tVar) = tProcRank;
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                    &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1 && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 0 ) <UINT_MAX ) // Bottom
            {
                tProcShare(tVar) = aProcNeighbor( 0 );
                tVar++;
            }
                        if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 )
                    &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1 && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 2 ) <UINT_MAX ) // Top
            {
                tProcShare(tVar) = aProcNeighbor( 2 );
                tVar++;
            }
                       if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1  &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1 && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 3 ) <UINT_MAX ) // Left
            {
                tProcShare(tVar) = aProcNeighbor( 3 );
                tVar++;
            }
                                   if ( tBasisPosition( 0 ) == tProcRange( 1 )  &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) - 1 &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1 && tBasisPosition( 2 ) <= tProcRange( 5 ) && aProcNeighbor( 1 ) <UINT_MAX ) // Right
            {
                tProcShare(tVar) = aProcNeighbor( 1 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbor( 4 ) <UINT_MAX ) // Back
            {
                tProcShare(tVar) = aProcNeighbor( 4 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1
                    && tBasisPosition( 1 ) <= tProcRange( 3 ) &&  tBasisPosition( 2 ) == tProcRange( 5 ) && aProcNeighbor( 5 ) <UINT_MAX ) // Front
            {
                tProcShare(tVar) = aProcNeighbor( 5 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1 && aProcNeighbor( 6 ) <UINT_MAX ) // Back bottom left
            {
                tProcShare(tVar) = aProcNeighbor( 6 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 9 ) <UINT_MAX ) // Back bottom
            {
                tProcShare(tVar) = aProcNeighbor( 9 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 10 ) <UINT_MAX ) // Back bottom right
            {
                tProcShare(tVar) = aProcNeighbor( 10 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1 &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 7 ) <UINT_MAX ) // Back left
            {
                tProcShare(tVar) = aProcNeighbor( 7 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1 &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 11 ) <UINT_MAX ) // Back right
            {
                tProcShare(tVar) = aProcNeighbor( 11 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 12 ) <UINT_MAX ) // Back top left
            {
                tProcShare(tVar) = aProcNeighbor( 12 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 13 ) <UINT_MAX ) // Back top
            {
                tProcShare(tVar) = aProcNeighbor( 13 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 4 ) - 1  && aProcNeighbor( 14 ) <UINT_MAX ) // Back top right
            {
                tProcShare(tVar) = aProcNeighbor( 14 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 8 ) <UINT_MAX ) // Bottom left
            {
                tProcShare(tVar) = aProcNeighbor( 8 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 15 ) <UINT_MAX ) // Bottom right
            {
                tProcShare(tVar) = aProcNeighbor( 15 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 16 ) <UINT_MAX ) // Top left
            {
                tProcShare(tVar) = aProcNeighbor( 16 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) >= tProcRange( 4 ) - 1  &&  tBasisPosition( 2 ) >= tProcRange( 5 )  && aProcNeighbor( 17 ) <UINT_MAX ) // Top right
            {
                tProcShare(tVar) = aProcNeighbor( 17 );
                tVar++;
            }
                        if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 18 ) <UINT_MAX ) // Front bottom left
            {
                tProcShare(tVar) = aProcNeighbor( 18 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 19 ) <UINT_MAX ) // Front bottom
            {
                tProcShare(tVar) = aProcNeighbor( 19 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 2 ) - 1
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 20 ) <UINT_MAX ) // Front bottom right
            {
                tProcShare(tVar) = aProcNeighbor( 20 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1 &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 21 ) <UINT_MAX ) // Front left
            {
                tProcShare(tVar) = aProcNeighbor( 21 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) >= tProcRange( 2 ) - 1 &&  tBasisPosition( 1 ) <= tProcRange( 3 )
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 22 ) <UINT_MAX ) // Front right
            {
                tProcShare(tVar) = aProcNeighbor( 22 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 0 ) - 1 &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 23 ) <UINT_MAX ) // Front top left
            {
                tProcShare(tVar) = aProcNeighbor( 23 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) >= tProcRange( 0 ) - 1 && tBasisPosition( 0 ) <= tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 24 ) <UINT_MAX ) // Front top
            {
                tProcShare(tVar) = aProcNeighbor( 24 );
                tVar++;
            }
            if ( tBasisPosition( 0 ) == tProcRange( 1 ) &&  tBasisPosition( 1 ) == tProcRange( 3 ) 
                     &&  tBasisPosition( 2 ) == tProcRange( 5 )  && aProcNeighbor( 25 ) <UINT_MAX ) // Front top right
            {
                tProcShare(tVar) = aProcNeighbor( 25 );
                tVar++;
            }
        }
    }
    tProcShare.resize(tVar,1);
    return tProcShare;
}

Mat<uint>
Lagrange_Basis::give_parents_of_basis(
        Mat<uint> const & aBasisList,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    //Level of basis
    uint tLevelOfBasis = 0;
    //"Parent" of basis
    uint tPossibleParentOfBasis = UINT_MAX;
    //Copy first the possible basis and then check if one of them has a "parent"
    Mat<uint> tParentOfBasis = aBasisList;
    //Determine the lowest possible basis id (Check for "parents" )
    Mat<uint> tLowestBasis;
    //Loop over all possible basis
    for ( uint i = 0; i < tParentOfBasis.length(); i++ )
    {
        tLevelOfBasis = give_basis_level( tParentOfBasis( i ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        tPossibleParentOfBasis = tParentOfBasis( i );
        while( tLevelOfBasis >  0 && tPossibleParentOfBasis < UINT_MAX )
        {
            //Parent of Basis is there
            tParentOfBasis( i ) = tPossibleParentOfBasis;
            //Check for next parent
            tPossibleParentOfBasis = give_basis_of_parent( tPossibleParentOfBasis, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        }
    }
    return tParentOfBasis;
}
