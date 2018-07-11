/*
 * cl_Lagrange_Element.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#include "cl_Lagrange_Element.hpp"
using namespace moris;

Mat<uint>
Lagrange_Element::give_basis_of_element(
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Basis functions of an element
    Mat<uint> tBasis( pow( aPolynomial + 1 , aModelDim ) , 1 );
    // Computes the level of the element
    uint tElementLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Determines the position of the element
    Mat<uint> tElementPosition = mBaseElement.give_position_of_element( aElementId, aModelDim, aNumberOfElementsPerDirection );
    uint tNumberOfBasis = 0;
    //Temporary variable for the for loop
    uint tVar = 0;
    if ( tElementLevel != 0 )
    {
        //Count basis functions until tBasis_level-1
        tNumberOfBasis = Lagrange_Basis::give_number_of_basis( tElementLevel - 1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    }
    // Temporary variable for simple coding style (Computing of number of basis in x-direction)
    uint tPartialPosition = pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 ) * aPolynomial + 1;
    if ( aModelDim == 1 )
    {
        if ( aPolynomial == 1 )
        {
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial     + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1 + tNumberOfBasis;
        }
        else if ( aPolynomial == 2 )
        {
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial     + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1 + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) * aPolynomial + 2 + tNumberOfBasis;
        }
        else if ( aPolynomial == 3 )
        {
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial     + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1 + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) * aPolynomial + 2 + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) * aPolynomial + 3 + tNumberOfBasis;
        }
        else
        {
            // The basis functions will be calculated, if higher order polynomials are needed
            for (uint i = 0; i < aPolynomial + 1; i++ )
            {
                tBasis(tVar) = tElementPosition( 0 ) * aPolynomial + i + tNumberOfBasis;;
                tVar++;
            }
        }
    }
    else if ( aModelDim == 2 )
    {
        if ( aPolynomial == 1 )
        {
            // Basis functions can be determined with the position of the element
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial
                    + ( tElementPosition( 1 ) * aPolynomial )     * tPartialPosition + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1
                    + ( tElementPosition( 1 ) * aPolynomial )     * tPartialPosition + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) * aPolynomial
                    + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) * aPolynomial + 1
                    + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition + tNumberOfBasis;
        }
        else if (aPolynomial == 2)
        {
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial
                    + ( tElementPosition( 1 ) * aPolynomial )     * tPartialPosition + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1
                    + ( tElementPosition( 1 ) * aPolynomial )     * tPartialPosition + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) * aPolynomial + 2
                    + ( tElementPosition( 1 ) * aPolynomial )     * tPartialPosition + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) * aPolynomial
                    + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 ) * aPolynomial + 1
                    + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 ) * aPolynomial + 2
                    + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 ) * aPolynomial
                    + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 ) * aPolynomial + 1
                    + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 8 ) = tElementPosition( 0 ) * aPolynomial + 2
                    + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition + tNumberOfBasis;
        }
        else
        {
            // The basis functions will be calculated, if higher order polynomials are needed
            for (uint i = 0; i<aPolynomial+1; i++)
            {
                for (uint j = 0; j<aPolynomial+1; j++)
                {
                    tBasis(tVar) = ( tElementPosition( 0 ) * aPolynomial + j )
                                 + ( tElementPosition( 1 ) * aPolynomial + i )* tPartialPosition + tNumberOfBasis;
                    tVar++;
                }
            }
        }
    }
    else if (aModelDim == 3)
    {
        //Temporary pre-calculation for simpler displaying
        uint tTempPreCalculation = tPartialPosition * ( pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 1 ) * aPolynomial + 1 );
        if (aPolynomial == 1)
        {
            tBasis( 0 ) = tElementPosition( 0 ) * aPolynomial     +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) * aPolynomial + 1 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) * aPolynomial     + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 ) * aPolynomial     +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 ) * aPolynomial + 1 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 ) * aPolynomial     + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
        }
        else if (aPolynomial == 2)
        {
            tBasis( 0 )  = tElementPosition( 0 ) * aPolynomial +       tElementPosition( 1 ) * aPolynomial       * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 1 )  = tElementPosition( 0 ) * aPolynomial + 1 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 2 )  = tElementPosition( 0 ) * aPolynomial + 2 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 3 )  = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 4 )  = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 5 )  = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 6 )  = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 7 )  = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 8 )  = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  +   tElementPosition( 2 ) * aPolynomial       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 9 )  = tElementPosition( 0 ) * aPolynomial +       tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 10 ) = tElementPosition( 0 ) * aPolynomial + 1 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 11 ) = tElementPosition( 0 ) * aPolynomial + 2 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 12 ) = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 13 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 14 ) = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 15 ) = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 16 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 17 ) = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 18 ) = tElementPosition( 0 ) * aPolynomial +       tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 19 ) = tElementPosition( 0 ) * aPolynomial + 1 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 20 ) = tElementPosition( 0 ) * aPolynomial + 2 +   tElementPosition( 1 ) * aPolynomial       * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 21 ) = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 22 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 23 ) = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 24 ) = tElementPosition( 0 ) * aPolynomial +     ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 25 ) = tElementPosition( 0 ) * aPolynomial + 1 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 26 ) = tElementPosition( 0 ) * aPolynomial + 2 + ( tElementPosition( 1 ) * aPolynomial + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) * aPolynomial + 2 ) * tTempPreCalculation + tNumberOfBasis;
        }
        else
        {
            // The basis functions will be calculated, if higher order polynomials are needed
            for ( uint k = 0; k < aPolynomial + 1; k++ )
            {
                for ( uint i = 0; i < aPolynomial + 1; i++ )
                {
                    for ( uint j = 0; j < aPolynomial + 1; j++ )
                    {
                        tBasis(tVar) = ( tElementPosition( 0 ) * aPolynomial + j )
                                     + ( tElementPosition( 1 ) * aPolynomial + i )* tPartialPosition
                                     + ( tElementPosition( 2 ) * aPolynomial + k )* tTempPreCalculation + tNumberOfBasis;
                        tVar++;
                    }
                }
            }
        }
    }
    return tBasis;
}
