/*
 * cl_Hierachical_Mesh_Element.cpp
 *
 *  Created on: Dec 7, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Element.hpp"
using namespace moris;

Mat<uint>
Hierarchical_Mesh_Element::give_basis_of_element(
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
        tNumberOfBasis = Hierarchical_Mesh_Basis::give_number_of_basis( tElementLevel-1, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    }
    // Temporary variable for simple coding style (Computing of a part of the position of the element)
    uint tPartialPosition = pow( 2, tElementLevel )*aNumberOfElementsPerDirection( 0 ) + aPolynomial;
    if ( aModelDim == 1 )
    {
        if ( aPolynomial == 1 )
        {
            tBasis( 0 ) = tElementPosition( 0 ) ;
            tBasis( 1 ) = tElementPosition( 0 ) + 1;
        }
        else if ( aPolynomial == 2 )
        {
            tBasis( 0 ) = tElementPosition( 0 );
            tBasis( 1 ) = tElementPosition( 0 ) + 1;
            tBasis( 2 ) = tElementPosition( 0 ) + 2;
        }
        else if ( aPolynomial == 3 )
        {
            tBasis( 0 ) = tElementPosition( 0 );
            tBasis( 1 ) = tElementPosition( 0 ) + 1;
            tBasis( 2 ) = tElementPosition( 0 ) + 2;
            tBasis( 3 ) = tElementPosition( 0 ) + 3;
        }
        else
        {
            // The basis functions will be calculated, if higher order polynomials are needed
            for (uint i = 0; i < aPolynomial + 1; i++ )
            {
                tBasis(tVar) = ( tElementPosition( 0 ) + i );
                tVar++;
            }
        }
    }
    else if ( aModelDim == 2 )
    {
        if ( aPolynomial == 1 )
        {
            // Basis functions can be determined with the position of the element
            tBasis( 0 ) = tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
        }
        else if (aPolynomial == 2)
        {
            tBasis( 0 ) = tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) + 2 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 8 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
        }
        else if (aPolynomial == 3)
        {
            tBasis( 0 ) = tElementPosition( 0 )  +       tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 )  + 1 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 )  + 2 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 )  + 3 +   tElementPosition( 1 )       * tPartialPosition + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 )  +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 )  + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 )  + 2 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 )  + 3 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 8 ) = tElementPosition( 0 )  +     ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 9 ) =  tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 10 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 11 ) = tElementPosition( 0 ) + 3 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 12 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 3 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 13 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 3 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 14 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 3 ) * tPartialPosition + tNumberOfBasis;
            tBasis( 15 ) = tElementPosition( 0 ) + 3 + ( tElementPosition( 1 ) + 3 ) * tPartialPosition + tNumberOfBasis;
        }
        else
        {
            // The basis functions will be calculated, if higher order polynomials are needed
            for (uint i = 0; i<aPolynomial+1; i++)
            {
                for (uint j = 0; j<aPolynomial+1; j++)
                {
                    tBasis(tVar) = ( tElementPosition( 0 ) + j )
                                                                                    + ( tElementPosition( 1 ) + i )* tPartialPosition + tNumberOfBasis;
                    tVar++;
                }
            }
        }
    }
    else if (aModelDim == 3)
    {
        //Temporary pre-calculation for simpler displaying
        uint tTempPreCalculation = tPartialPosition *(pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 1 ) + aPolynomial );
        if (aPolynomial == 1)
        {
            tBasis( 0 ) = tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
        }
        else if (aPolynomial == 2)
        {
            tBasis( 0 ) = tElementPosition( 0 ) +        tElementPosition( 1 )       * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 1 ) = tElementPosition( 0 ) + 1 +    tElementPosition( 1 )       * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 2 ) = tElementPosition( 0 ) + 2 +    tElementPosition( 1 )       * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 3 ) = tElementPosition( 0 ) +      ( tElementPosition( 1 ) + 1 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 4 ) = tElementPosition( 0 ) + 1 +  ( tElementPosition( 1 ) + 1 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 5 ) = tElementPosition( 0 ) + 2 +  ( tElementPosition( 1 ) + 1 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 6 ) = tElementPosition( 0 ) +      ( tElementPosition( 1 ) + 2 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 7 ) = tElementPosition( 0 ) + 1 +  ( tElementPosition( 1 ) + 2 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 8 ) = tElementPosition( 0 ) + 2 +  ( tElementPosition( 1 ) + 2 ) * tPartialPosition  +   tElementPosition( 2 )       * tTempPreCalculation + tNumberOfBasis;
            tBasis( 9 ) =  tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 10 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 11 ) = tElementPosition( 0 ) + 2 +   tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 12 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 13 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 14 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 15 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 16 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 17 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 1 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 18 ) = tElementPosition( 0 ) +       tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) +2 )  * tTempPreCalculation + tNumberOfBasis;
            tBasis( 19 ) = tElementPosition( 0 ) + 1 +   tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) +2 )  * tTempPreCalculation + tNumberOfBasis;
            tBasis( 20 ) = tElementPosition( 0 ) + 2 +   tElementPosition( 1 )       * tPartialPosition  + ( tElementPosition( 2 ) +2 )  * tTempPreCalculation + tNumberOfBasis;
            tBasis( 21 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 22 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 23 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 1 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 24 ) = tElementPosition( 0 ) +     ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 25 ) = tElementPosition( 0 ) + 1 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
            tBasis( 26 ) = tElementPosition( 0 ) + 2 + ( tElementPosition( 1 ) + 2 ) * tPartialPosition  + ( tElementPosition( 2 ) + 2 ) * tTempPreCalculation + tNumberOfBasis;
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
                        tBasis(tVar) = ( tElementPosition( 0 ) + j )
                                                                                        + ( tElementPosition( 1 ) + i )* tPartialPosition
                                                                                        + ( tElementPosition( 2 ) + k )* tTempPreCalculation + tNumberOfBasis;
                        tVar++;
                    }
                }
            }
        }
    }
    return tBasis;
}

//--------------------------------------------------------------------------------

/* @TODO give_middlecoordinate_from_element: simplify formula */
Mat<real>
Hierarchical_Mesh_Element::give_middlecoordinate_from_element(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<real> const & aModelDimensions,
        Mat<real> const & aModelDimensions_Offset) const
{
    //Coordinates of the middle coordinate of the element
    Mat<real> tCoordinates( aModelDim , 1 , 0 );
    // Compute the level of the element
    uint tElementLevel = mBaseElement.give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Get the Position of the element
    Mat<uint> tElementPosition = mBaseElement.give_position_of_element( aElementId, aModelDim, aNumberOfElementsPerDirection );
    //Calculate the Coordinates in the middle  of the element with respect to the dimensions and the offset
    for ( uint i = 0; i < aModelDim; i++ )
    {
        //XYZ-coordinates + offset (if needed);
        tCoordinates( i ) = ( tElementPosition( i ) * aModelDimensions( i ) ) / ( aNumberOfElementsPerDirection( i ) * pow( 2, tElementLevel ) )
                                        + aModelDimensions_Offset( i ) + aModelDimensions( i ) / (aNumberOfElementsPerDirection( i ) * pow( 2, tElementLevel ) * 2 );
    }
    return tCoordinates;
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Element::give_element_for_coordinate_on_level_zero(
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<real> const & aDomainRange,
        Mat<real> const & aDomainOffset,
        Mat<real> const & aPointOfOrigin,
        Mat<real> const & aCoordinate) const
{
    Mat<uint> tIJKPosition( aModelDim , 1 , 0 );
    uint tElement = UINT_MAX;
    if ( aModelDim == 1 )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = (aCoordinate( 0 ) - aDomainOffset( 0 ) ) / aDomainRange( 0 );
        // Calculates the with the ratio the position in the domain
        tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) );
        // Defines the element for this position
        tElement = mBaseElement.give_element_of_position( 0, aModelDim, aNumberOfElementsPerDirection, tIJKPosition );
    }
    else if ( aModelDim == 2 )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = (aCoordinate( 0 )-aDomainOffset( 0 ))/aDomainRange( 0 );
        // Calculates the ratio within the domain in y-direction
        real tYRatio = (aCoordinate( 1 )-aDomainOffset( 1 ))/aDomainRange( 1 );
        // Calculates the with the ratio the position in the domain
        tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) );
        tIJKPosition( 1 ) = floor( tYRatio * aNumberOfElementsPerDirection( 1 ) );
        // Defines the element for this position
        tElement = mBaseElement.give_element_of_position( 0, aModelDim, aNumberOfElementsPerDirection, tIJKPosition );
    }
    else if ( aModelDim == 3 )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = (aCoordinate( 0 )-aDomainOffset( 0 ))/aDomainRange( 0 );
        // Calculates the ratio within the domain in y-direction
        real tYRatio = (aCoordinate( 1 )-aDomainOffset( 1 ))/aDomainRange( 1 );
        // Calculates the ratio within the domain in z-direction
        real tZRatio = (aCoordinate( 2 )-aDomainOffset( 2 ))/aDomainRange( 2 );
        // Calculates the with the ratio the position in the domain
        tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) );
        tIJKPosition( 1 ) = floor( tYRatio * aNumberOfElementsPerDirection( 1 ) );
        tIJKPosition( 2 ) = floor( tZRatio * aNumberOfElementsPerDirection( 2 ) );
        // Defines the element for this position
        tElement = mBaseElement.give_element_of_position( 0, aModelDim, aNumberOfElementsPerDirection, tIJKPosition );
    }
    return tElement;
}

//--------------------------------------------------------------------------------

uint
Hierarchical_Mesh_Element::give_active_element_for_coordinate(
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<real> const & aDomainRange,
        Mat<real> const & aDomainOffset,
        Mat<real> const & aPointOfOrigin,
        Mat<real> const & aCoordinate,
        BoostBitset const & aElementActive) const
{
    uint  tElementId = UINT_MAX;
    uint tLevel = 0;
    //Temporary variable to determine the element id
    uint tElementDummy = 0;
    //Position of element with i,j,k coordinates
    Mat<uint> tIJKPosition( aModelDim , 1 , 0 );
    bool tInDomainCheck = false; // Flag is used, to check if the coordinate is in the domain (Switch for MORIS_ASSERT)
    if ( aModelDim == 1 &&  aCoordinate( 0 ) > aPointOfOrigin( 0 )
            && aCoordinate( 0 ) < ( aDomainRange( 0 ) + 2 * aDomainOffset( 0 ) + aPointOfOrigin( 0 ) ) )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = (aCoordinate( 0 )-aDomainOffset( 0 ))/aDomainRange( 0 );
        while ( tElementId == UINT_MAX )
        {
            // Calculates the with the ratio the position in the domain
            tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) * pow( 2, tLevel ) );
            // Defines the element for this position
            tElementDummy = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tIJKPosition );
            // If the element is active, the right element is found
            if ( tElementDummy < aElementActive.size() && aElementActive.test( tElementDummy ) == 1 )
            {
                tElementId = tElementDummy;
                tInDomainCheck = true;
            }
            else if ( tElementDummy > aElementActive.size() )
            {
                break;
            }
            tLevel++;
        }
    }
    else if ( aModelDim == 2 &&  aCoordinate( 0 ) > aPointOfOrigin( 0 )
            && aCoordinate( 0 ) < ( aDomainRange( 0 ) + 2 * aDomainOffset( 0 ) + aPointOfOrigin( 0 ) )
            && aCoordinate( 1 ) > aPointOfOrigin( 1 ) && aCoordinate( 1 ) < ( aDomainRange( 1 ) + 2 * aDomainOffset( 1 ) + aPointOfOrigin( 1 ) ) )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = ( aCoordinate( 0 ) - aDomainOffset( 0 ) ) / aDomainRange( 0 );
        // Calculates the ratio within the domain in y-direction
        real tYRatio = ( aCoordinate( 1 ) - aDomainOffset( 1 ) ) / aDomainRange( 1 );
        while ( tElementId == UINT_MAX )
        {
            // Calculates the with the ratio the position in the domain
            tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) * pow( 2, tLevel ) );
            tIJKPosition( 1 ) = floor( tYRatio * aNumberOfElementsPerDirection( 1 ) * pow( 2, tLevel ) );
            // Defines the element for this position
            tElementDummy = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tIJKPosition );
            // If the element is active, the right element is found
            if ( tElementDummy < aElementActive.size() && aElementActive.test( tElementDummy ) == 1 )
            {
                tElementId = tElementDummy;
                tInDomainCheck = true;
            }
            else if ( tElementDummy > aElementActive.size() )
            {
                break;
            }
            tLevel++;
        }
    }
    else if ( aModelDim == 3 && aCoordinate( 0 ) > aPointOfOrigin( 0 )
            && aCoordinate( 0 ) < ( aDomainRange( 0 ) + 2 * aDomainOffset( 0 ) + aPointOfOrigin( 0 ) )
            && aCoordinate( 1 ) > aPointOfOrigin( 1 ) && aCoordinate( 1 ) < ( aDomainRange( 1 ) + 2 * aDomainOffset( 1 ) + aPointOfOrigin( 1 ) )
            && aCoordinate( 2 ) > aPointOfOrigin( 2 ) && aCoordinate( 2 ) < ( aDomainRange( 2 ) + 2 * aDomainOffset( 2 ) + aPointOfOrigin( 2 ) ) )
    {
        // Calculates the ratio within the domain in x-direction
        real tXRatio = (aCoordinate( 0 )-aDomainOffset( 0 ))/aDomainRange( 0 );
        // Calculates the ratio within the domain in y-direction
        real tYRatio = (aCoordinate( 1 )-aDomainOffset( 1 ))/aDomainRange( 1 );
        // Calculates the ratio within the domain in z-direction
        real tZRatio = (aCoordinate( 2 )-aDomainOffset( 2 ))/aDomainRange( 2 );
        while ( tElementId == UINT_MAX )
        {
            // Calculates the with the ratio the position in the domain
            tIJKPosition( 0 ) = floor( tXRatio * aNumberOfElementsPerDirection( 0 ) * pow( 2, tLevel ) );
            tIJKPosition( 1 ) = floor( tYRatio * aNumberOfElementsPerDirection( 1 ) * pow( 2, tLevel ) );
            tIJKPosition( 2 ) = floor( tZRatio * aNumberOfElementsPerDirection( 2 ) * pow( 2, tLevel ) );
            // Defines the element for this position
            tElementDummy = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tIJKPosition);
            // If the element is active, the right element is found
            if ( tElementDummy < aElementActive.size() && aElementActive.test( tElementDummy ) == 1)
            {
                tElementId = tElementDummy;
                tInDomainCheck = true;
            }
            else if ( tElementDummy > aElementActive.size() )
            {
                break;
            }
            tLevel++;
        }
    }
    MORIS_ERROR(tInDomainCheck == true,"Coordinate is not in the domain");
    return tElementId;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Element::give_active_face_neighbor_of_element(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aLevel,
        BoostBitset const & aElementActive) const
{
    // Neighbors of the element
    Mat<uint> tNeighbor = mBaseElement.give_neighbor_of_element( aElementId, aModelDim, 1, aNumberOfElementsPerDirection );
    Mat<uint> tPossibleNeighbors;
    // Estimation of the possible number of neighbors
    //( First col is element neighbor id, second col is Ordinal of aElementId to Neighbor and third col is ordinal from neighbor to aElementId)
    Mat<uint> tActiveNeighbors( 2 * aModelDim * pow( pow( 2, aLevel ) , aModelDim ) , 3 , UINT_MAX );
    //Temporary variables for loop

    // counter for tActiveNeighbors
    uint tVar = 0;

    // counter for tPossibleNeighbors (reading)
    uint tVara = 0;

    // counter for tPossibleNeighbors (writing)
    uint tVarb = 1;
    uint tParent = 0;
    uint tLevel = 0;
    Mat<uint> tChildren;
    Mat<uint> tNeighbors;
    Mat<uint> tNeighborsOrdinal;
    Mat<uint> tOrdinalOfNeighbors;
    Mat<uint> tChildNeighbors;
    Mat<uint> tChildNeighborsOrdinal;
    Mat<uint> tOrdinalOfChildNeighbors;
    uint tSwitch = 0;
    if ( aModelDim == 1 )
    {
        // Neighbors in 1D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors1D = { { 0, 2 } };
        tNeighbors = tNeighbors1D;
        // Neighbors in 1D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal1D = { { 3, 1 } };
        tNeighborsOrdinal = tNeighborsOrdinal1D;
        // Neighbors ordinal to the element Id in 1D
        Mat<uint> tOrdinalOfNeighbors1D = { { 1, 3 } };
        tOrdinalOfNeighbors = tOrdinalOfNeighbors1D;
        // Possible childrens for each neighbor
        Mat<uint> tChildNeighbors1D = { { 1 },{ 0 } };
        tChildNeighbors = tChildNeighbors1D;
        // Ordinal to childrens for each neighbor
        Mat<uint> tChildNeighborsOrdinal1D = { { 3 }, { 1 } };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal1D;
        // Ordinal from childrens for each neighbor
        Mat<uint> tOrdinalOfChildNeighbor1D = { { 1 }, { 3 } };
        tOrdinalOfChildNeighbors = tOrdinalOfChildNeighbor1D;
    }
    else if ( aModelDim == 2 )
    {
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors2D = { {1, 3, 5, 7 } };
        tNeighbors = tNeighbors2D;
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal2D = { {0, 3, 1, 2 } };
        tNeighborsOrdinal = tNeighborsOrdinal2D;
        // Neighbors ordinal to the element Id in 1D
        Mat<uint> tOrdinalOfNeighbors2D = { {2, 1, 3, 0 } };
        tOrdinalOfNeighbors = tOrdinalOfNeighbors2D;
        // Possible childrens for each neighbor
        Mat<uint> tChildNeighbors2D = { {2, 3}, {1, 3}, {0, 2}, {0, 1} };
        tChildNeighbors = tChildNeighbors2D;
        // Ordinal to childrens for each neighbor
        Mat<uint> tChildNeighborsOrdinal2D = { {0, 0}, {3, 3}, {1, 1}, {2, 2} };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal2D;
        // Ordinal from childrens for each neighbor
        Mat<uint> tOrdinalOfChildNeighbor2D = { {2, 2}, {1, 1}, {3, 3}, {0, 0} };
        tOrdinalOfChildNeighbors = tOrdinalOfChildNeighbor2D;
    }
    else if ( aModelDim == 3 )
    {
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors3D = { {4, 10, 12, 14, 16, 22} };
        tNeighbors = tNeighbors3D;
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal3D = { {4, 0, 3, 1, 2, 5} };
        tNeighborsOrdinal = tNeighborsOrdinal3D;
        // Neighbors ordinal to the element Id in 1D
        Mat<uint> tOrdinalOfNeighbors3D = { {5, 2, 1, 3, 0, 4} };
        tOrdinalOfNeighbors = tOrdinalOfNeighbors3D;
        // Possible childrens for each neighbor
        Mat<uint> tChildNeighbors3D = { {4, 5, 6, 7}, {2, 3, 6, 7}, {1, 3, 5, 7},
                {0, 2, 4, 6}, {0, 1, 4, 5}, {0, 1, 2, 3} };
        tChildNeighbors = tChildNeighbors3D;
        // Ordinal to childrens for each neighbor
        Mat<uint> tChildNeighborsOrdinal3D = { {4, 4, 4, 4}, {0, 0, 0, 0}, {3, 3, 3, 3},
                {1, 1, 1, 1}, {2, 2, 2, 2}, {5, 5, 5, 5} };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal3D;
        // Ordinal from childrens for each neighbor
        Mat<uint> tOrdinalOfChildNeighbor3D = { {5, 5, 5, 5}, {2, 2, 2, 2}, {1, 1, 1, 1},
                {3, 3, 3, 3}, {0, 0, 0, 0}, {4, 4, 4, 4} };
        tOrdinalOfChildNeighbors = tOrdinalOfChildNeighbor3D;
    }
    for ( uint i = 0; i < tNeighbors.length(); i++ )
    {
        MORIS_ASSERT( tNeighbor(tNeighbors( i ) ) < aElementActive.size(), "UINT_MAX cannot be an neighbor element!");
        //It is a neighbor, if active
        if ( aElementActive.test( tNeighbor(tNeighbors( i ) ) ) == 1 )
        {
            tActiveNeighbors( tVar, 0 ) = tNeighbor( tNeighbors( i ) );
            tActiveNeighbors( tVar, 1 ) = tNeighborsOrdinal( i );
            tActiveNeighbors( tVar, 2 ) = tOrdinalOfNeighbors( i );
            tVar++;
        }
        else
        {
            // If not active, then check parent and children
            tLevel = mBaseElement.give_element_level( tNeighbor( tNeighbors( i ) ), aModelDim, aNumberOfElementsPerDirection );
            tPossibleNeighbors.set_size( 2 * aModelDim * pow(2,aLevel) , 1 , 0 );
            tPossibleNeighbors( 0 ) =  tNeighbor(tNeighbors( i ));
            tVara = 0;
            tVarb = 1;
            tSwitch = 0;
            while( tPossibleNeighbors( tVara ) > 0 )
            {
                tParent = mBaseElement.give_parent_of_element( tPossibleNeighbors( tVara ), aModelDim, aNumberOfElementsPerDirection );
                if ( tParent < UINT_MAX && aElementActive.test(tParent) == 1 )
                {
                    tActiveNeighbors( tVar, 0 ) = tParent;
                    tActiveNeighbors( tVar, 1 ) = tNeighborsOrdinal( i );
                    tActiveNeighbors( tVar, 2 ) = tOrdinalOfNeighbors( i );
                    tVar++;
                    // If parent is found, no need to check childrens
                    tSwitch = 1;
                }
                else if (  tParent < UINT_MAX )
                {
                    tPossibleNeighbors( tVarb ) =  tParent;
                    tVarb++;
                }
                tVara++;
            }
            tPossibleNeighbors.set_size( 2 * aModelDim * pow( 2, aLevel ) , 1 , 0 );
            tPossibleNeighbors( 0 ) = tNeighbor( tNeighbors( i ) );
            tVara = 0;
            tVarb = 1;
            if ( tLevel + 1 <= aLevel )
            {
                while( tPossibleNeighbors( tVara ) > 0 && tSwitch == 0 )
                {
                    // Give the children of the parent element
                    tChildren = mBaseElement.give_children_of_element( tPossibleNeighbors( tVara ), aModelDim, aNumberOfElementsPerDirection );
                    tLevel = mBaseElement.give_element_level( tChildren( 0 ), aModelDim, aNumberOfElementsPerDirection );
                    if ( aElementActive.test( tChildren( tChildNeighbors( i, 0 ) ) ) == 1 )
                    {
                        tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 0) );
                        tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 0 );
                        tActiveNeighbors( tVar, 2 ) = tOrdinalOfChildNeighbors( i, 0 );
                        tVar++;
                    }
                    else if ( tLevel < aLevel )
                    {
                        tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 0 ) );
                        tVarb++;
                    }
                    if ( aElementActive.test( tChildren( tChildNeighbors( i, 1 ) ) ) == 1 )
                    {
                        tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 1 ) );
                        tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 1 );
                        tActiveNeighbors( tVar, 2 ) = tOrdinalOfChildNeighbors( i, 1 );
                        tVar++;
                    }
                    else if ( tLevel < aLevel )
                    {
                        tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 1 ) );
                        tVarb++;
                    }
                    if ( aModelDim == 3)
                    {
                        if ( aElementActive.test( tChildren( tChildNeighbors( i, 2) ) ) == 1 )
                        {
                            tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 2 ) );
                            tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 2 );
                            tActiveNeighbors( tVar, 2 ) = tOrdinalOfChildNeighbors( i, 2 );
                            tVar++;
                        }
                        else if ( tLevel < aLevel )
                        {
                            tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 2 ) );
                            tVarb++;
                        }
                        if ( aElementActive.test( tChildren( tChildNeighbors( i, 3 ) ) ) == 1 )
                        {
                            tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 3 ) );
                            tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 3 );
                            tActiveNeighbors( tVar, 2 ) = tOrdinalOfChildNeighbors( i, 3 );
                            tVar++;
                        }
                        else if ( tLevel < aLevel )
                        {
                            tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 3 ) );
                            tVarb++;
                        }
                    }
                    tVara++;
                }
            }
        }
    }
    tActiveNeighbors.resize( tVar , 3 );
    return tActiveNeighbors;
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Element::give_active_neighbor_of_element(
        uint const & aElement,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aLevel,
        BoostBitset const & aElementActive) const
{
    // Neighbors of the element
    Mat<uint> tNeighbor = mBaseElement.give_neighbor_of_element( aElement, aModelDim, 1, aNumberOfElementsPerDirection );
//    tNeighbor.print("tNeighbor");
    Mat<uint> tPossibleNeighbors;
    // Estimation of the possible number of neighbors
    Mat<uint> tActiveNeighbors( 2 * aModelDim * pow( pow( 2, aLevel ) , aModelDim ) , 2 , UINT_MAX );
    //Temporary variables for loop
    uint tVar = 0;
    uint tVara = 0;
    uint tVarb = 1;
    uint tParent = 0;
    uint tLevel = 0;
    Mat<uint> tChildren;
    Mat<uint> tNeighbors;
    Mat<uint> tNeighborsOrdinal;
    Mat<uint> tChildNeighbors;
    Mat<uint> tChildNeighborsOrdinal;
    uint tSwitch = 0;
    if ( aModelDim == 1 )
    {
        // Neighbors in 1D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors1D = { { 0, 2 } };
        tNeighbors = tNeighbors1D;
        // Neighbors in 1D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal1D = { { 3, 1 } };
        tNeighborsOrdinal = tNeighborsOrdinal1D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighbors1D = { { 1 },{ 0 } };
        tChildNeighbors = tChildNeighbors1D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighborsOrdinal1D = { { 3 }, { 0 } };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal1D;
    }
    else if ( aModelDim == 2 )
    {
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors2D = { {0, 1, 2, 3, 5, 6, 7, 8 } };
        tNeighbors = tNeighbors2D;
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal2D = { {0, 1, 2, 3, 5, 6, 7, 8 } };
        tNeighborsOrdinal = tNeighborsOrdinal2D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighbors2D = { {3, UINT_MAX }, {2, 3}, {2, UINT_MAX },
                {1, 3}, {0, 2}, {1, UINT_MAX }, {0, 1}, {0, UINT_MAX } };
        tChildNeighbors = tChildNeighbors2D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighborsOrdinal2D = { {0, UINT_MAX }, {1, 1}, {2, UINT_MAX },
                {3, 3}, {5, 5}, {6, UINT_MAX }, {7, 7}, {8, UINT_MAX }, };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal2D;
    }
    else if ( aModelDim == 3 )
    {
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighbors3D = { {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14,
                15, 16, 17, 18, 19, 20,21, 22, 23, 24, 25, 26} };
        tNeighbors = tNeighbors3D;
        // Neighbors in 2D (positions extracted from give_neighbor_of_element)
        Mat<uint> tNeighborsOrdinal3D = { {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                12, 14, 15, 16, 17, 18, 19, 20,21, 22, 23, 24, 25, 26} };
        tNeighborsOrdinal = tNeighborsOrdinal3D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighbors3D = { {7, UINT_MAX, UINT_MAX, UINT_MAX}, {6, 7, UINT_MAX, UINT_MAX},
                {6, UINT_MAX, UINT_MAX, UINT_MAX}, {5, 7, UINT_MAX, UINT_MAX}, {4, 5, 6, 7},
                {4, 6, UINT_MAX, UINT_MAX}, {5, UINT_MAX, UINT_MAX, UINT_MAX}, {4, 5, UINT_MAX, UINT_MAX},
                {4, UINT_MAX, UINT_MAX, UINT_MAX}, {3, 7, UINT_MAX, UINT_MAX}, {2, 3, 6, 7},
                {2, 6, UINT_MAX, UINT_MAX}, {1, 3, 5, 7}, {0, 2, 4, 6}, {1, 5, UINT_MAX, UINT_MAX},
                {0, 1, 4, 5}, {0, 4, UINT_MAX, UINT_MAX}, {3, UINT_MAX, UINT_MAX, UINT_MAX},
                {2, 3, UINT_MAX, UINT_MAX}, {2, UINT_MAX, UINT_MAX, UINT_MAX}, {1, 3, UINT_MAX, UINT_MAX},
                {0, 1, 2, 3}, {0, 2, UINT_MAX, UINT_MAX}, {1, UINT_MAX, UINT_MAX, UINT_MAX},
                {0, 1, UINT_MAX, UINT_MAX}, {0, UINT_MAX, UINT_MAX, UINT_MAX} };
        tChildNeighbors = tChildNeighbors3D;
        // index for possible children for each neighbor
        Mat<uint> tChildNeighborsOrdinal3D = { {0, UINT_MAX, UINT_MAX, UINT_MAX}, {1, 1, UINT_MAX, UINT_MAX},
                {2, UINT_MAX, UINT_MAX, UINT_MAX}, {3, 3, UINT_MAX, UINT_MAX}, {4, 4, 4, 4},
                {5, 5, UINT_MAX, UINT_MAX}, {6, UINT_MAX, UINT_MAX, UINT_MAX}, {7, 7, UINT_MAX, UINT_MAX},
                {8, UINT_MAX, UINT_MAX, UINT_MAX}, {9, 9, UINT_MAX, UINT_MAX}, {10, 10, 10, 10},
                {11, 11, UINT_MAX, UINT_MAX}, {12, 12, 12, 12}, {14, 14, 14, 14},
                {15, 15, UINT_MAX, UINT_MAX}, {16, 16, 16, 16}, {17, 17, UINT_MAX, UINT_MAX},
                {18, UINT_MAX, UINT_MAX, UINT_MAX}, {19, 19, UINT_MAX, UINT_MAX}, {20, UINT_MAX, UINT_MAX, UINT_MAX},
                {21, 21, UINT_MAX, UINT_MAX}, {22, 22, 22, 22}, {23, 23, UINT_MAX, UINT_MAX},
                {24, UINT_MAX, UINT_MAX, UINT_MAX}, {25, 25, UINT_MAX, UINT_MAX}, {26, UINT_MAX, UINT_MAX, UINT_MAX} };
        tChildNeighborsOrdinal = tChildNeighborsOrdinal3D;
    }
    for ( uint i = 0; i < tNeighbors.length(); i++ )
    {
        MORIS_ASSERT( tNeighbor(tNeighbors( i ) ) < aElementActive.size(), "UINT_MAX cannot be an neighbor element!");
        //It is a neighbor, if active
        if ( aElementActive.test( tNeighbor(tNeighbors( i ) ) ) == 1 )
        {
            tActiveNeighbors( tVar, 0 ) = tNeighbor( tNeighbors( i ) );
            tActiveNeighbors( tVar, 1 ) = tNeighborsOrdinal( i );
            tVar++;
        }
        else
        {
            // If not active, then check parent and children
            tLevel = mBaseElement.give_element_level( tNeighbor( tNeighbors( i ) ), aModelDim, aNumberOfElementsPerDirection );
            tPossibleNeighbors.set_size( 2 * aModelDim * pow(2,aLevel) , 1 , 0 );
            tPossibleNeighbors( 0 ) =  tNeighbor(tNeighbors( i ));
            tVara = 0;
            tVarb = 1;
            tSwitch = 0;
            while( tPossibleNeighbors( tVara ) > 0 )
            {
//                tPossibleNeighbors.print("tPossibleNeighbors");
                tParent = mBaseElement.give_parent_of_element( tPossibleNeighbors( tVara ), aModelDim, aNumberOfElementsPerDirection );
                if ( tParent < UINT_MAX && aElementActive.test(tParent) == 1 )
                {
                    tActiveNeighbors( tVar, 0 ) = tParent;
                    tActiveNeighbors( tVar, 1 ) = tNeighborsOrdinal( i );
                    tVar++;
                    // If parent is found, no need to check childrens
                    tSwitch = 1;
                }
                else if (  tParent < UINT_MAX )
                {
                    tPossibleNeighbors( tVarb ) =  tParent;
                    tVarb++;
                }
                tVara++;
            }
            tPossibleNeighbors.set_size( 2 * aModelDim * pow( 2, aLevel ) , 1 , 0 );
            tPossibleNeighbors( 0 ) = tNeighbor( tNeighbors( i ) );
            tVara = 0;
            tVarb = 1;
            if ( tLevel + 1 <= aLevel )
            {
                while( tPossibleNeighbors( tVara ) > 0 && tSwitch == 0 )
                {
                    // Give the children of the parent element
                    tChildren = mBaseElement.give_children_of_element( tPossibleNeighbors( tVara ), aModelDim, aNumberOfElementsPerDirection );
//                    tChildren.print("tChildren");
                    tLevel = mBaseElement.give_element_level( tChildren( 0 ), aModelDim, aNumberOfElementsPerDirection );
                    if ( tChildNeighbors( i, 0 ) < UINT_MAX && aElementActive.test( tChildren( tChildNeighbors( i, 0 ) ) ) == 1 )
                    {
                        tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 0) );
                        tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 0 );
                        tVar++;
                    }
                    else if ( tChildNeighbors( i, 0 ) < UINT_MAX && tLevel < aLevel )
                    {
                        tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 0 ) );
                        tVarb++;
                    }
                    if ( tChildNeighbors( i, 1 ) < UINT_MAX && aElementActive.test( tChildren( tChildNeighbors( i, 1 ) ) ) == 1 )
                    {
                        tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 1 ) );
                        tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 1 );
                        tVar++;
                    }
                    else if ( tChildNeighbors( i, 1 ) < UINT_MAX && tLevel < aLevel )
                    {
                        tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 1 ) );
                        tVarb++;
                    }
                    if ( aModelDim == 3)
                    {
                        if ( tChildNeighbors( i, 2 ) < UINT_MAX && aElementActive.test( tChildren( tChildNeighbors( i, 2) ) ) == 1 )
                        {
                            tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 2 ) );
                            tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 2 );
                            tVar++;
                        }
                        else if ( tChildNeighbors( i, 2 ) < UINT_MAX && tLevel < aLevel )
                        {
                            tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 2 ) );
                            tVarb++;
                        }
                        if ( tChildNeighbors( i, 3 ) < UINT_MAX && aElementActive.test( tChildren( tChildNeighbors( i, 3 ) ) ) == 1 )
                        {
                            tActiveNeighbors( tVar, 0 ) = tChildren( tChildNeighbors( i, 3 ) );
                            tActiveNeighbors( tVar, 1 ) = tChildNeighborsOrdinal( i, 3 );
                            tVar++;
                        }
                        else if ( tChildNeighbors( i, 3 ) < UINT_MAX && tLevel < aLevel )
                        {
                            tPossibleNeighbors( tVarb ) =  tChildren( tChildNeighbors( i, 3 ) );
                            tVarb++;
                        }
                    }
                    tVara++;
                }
            }
        }
    }
    tActiveNeighbors.resize( tVar , 2 );
    return tActiveNeighbors;
}

//--------------------------------------------------------------------------------

/* @TODO replace pow(2, tElementLevel) with constant */
Mat<uint>
Hierarchical_Mesh_Element::give_neighbor_stencil_of_element(
        uint const & aFeatureResolution,
        uint const & aElementId,
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset & ElementActive)
{
    // Determines the position of the element
    Mat<uint> tElementPosition = mBaseElement.give_position_of_element(aElementId,aModelDim,aNumberOfElementsPerDirection);
    // Computes the element level
    uint tElementLevel = mBaseElement.give_element_level(aElementId,aModelDim,aNumberOfElementsPerDirection);
    // Position of a neighbor element;
    Mat<uint> tNeighborPosition( aModelDim , 1 , UINT_MAX );
    Mat<uint> tElementNeighbor;
    // Temporary variable for loop

    // counter for tElementNeighbor
    uint tVar = 0;
    uint tElement;
    int i = 0;
    int j = 0;
    int k = 0;
    MORIS_ASSERT( aModelDim > 1, "Stencil is not implemented for 1D");
    if ( aModelDim == 2 )
    {
        // Initialize with zeros
        tElementNeighbor.set_size(4,2*aFeatureResolution+3,0);
        tVar = 0;
        // Y-direction stays constant
        j = 0;
        // Loop over x-direction
        for (i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            // Neighbor needs to be in the range of the elements in each direction
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - aPolynomial)
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;   tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                // Computes the element id with the position
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1 )
                {
                    // Save the neighbor if he is active
                    tElementNeighbor(0,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // X-direction is constant
        i = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - aPolynomial)
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;   tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1 )
                {
                    tElementNeighbor(1,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            // Checks for elements in a diagonal direction from bottom left to top right
            i = j;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - aPolynomial)
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(2,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            // Checks for elements in a diagonal direction from top left to bottom right
            i = -j;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - aPolynomial)
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(3,tVar) = tElement;
                }
            }
            tVar++;
        }
    }
    else if ( aModelDim == 3 )
    {
        tElementNeighbor.set_size( 12, 2*aFeatureResolution+3 ,0);
        tVar = 0;
        // Z-direction is constant
        k = 0;
        // Y-direction is constant
        j = 0;
        // Loop over x-direction
        for ( i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - aPolynomial
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < aNumberOfElementsPerDirection( 2 ) * pow( 2, tElementLevel ) - aPolynomial)
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(0,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Z-direction is constant
        k = 0;
        // X-direction is constant
        i = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(1,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // y-direction is constant
        j = 0;
        // X-direction is constant
        i = 0;
        // Loop over z-direction
        for ( k = -((int)aFeatureResolution+1); k<=((int)aFeatureResolution+1); k++ )
        {
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(2,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Z-direction is constant
        k = 0;
        // Loop over x-direction
        for ( i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            // Checks elements in the x-,y-direction
            j = i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(3,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Z-direction is constant
        k = 0;
        // Loop over x-direction
        for ( i = -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            // Checks elements in the x-,y-direction
            j = -i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(4,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // X-direction is constant
        i = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            // Checks elements in the y-,z-direction
            k = j;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(5,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // X-direction is constant
        i = 0;
        // Loop over y-direction
        for ( j = -((int)aFeatureResolution+1); j<=((int)aFeatureResolution+1); j++ )
        {
            // Checks elements in the y-,z-direction
            k = -j;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(6,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // y-direction is constant
        j = 0;
        // Loop over x-direction
        for ( i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            // Checks elements in the x-,z-direction
            k = i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(7,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // y-direction is constant
        j = 0;
        // Loop over x-direction
        for ( i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            // Checks elements in the x-,z-direction
            k = -i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(8,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Loop over x-direction
        for ( i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++ )
        {
            // Loop for a diagonal in 3D
            k = i;
            j = i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(9,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Loop over x-direction
        for (i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            // Loop for a diagonal in 3D
            k = i;
            j = -i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(10,tVar) = tElement;
                }
            }
            tVar++;
        }
        tVar = 0;
        // Loop over x-direction
        for (i= -((int)aFeatureResolution+1); i<=((int)aFeatureResolution+1); i++)
        {
            // Loop for a diagonal in 3D
            k = -i;
            j = i;
            if ( tElementPosition( 0 ) + i >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 0 ) + i < (aNumberOfElementsPerDirection( 0 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 1 ) + j < (aNumberOfElementsPerDirection( 1 ) - aPolynomial) * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k >= aPolynomial * pow( 2, tElementLevel )
            && tElementPosition( 2 ) + k < (aNumberOfElementsPerDirection( 2 ) - aPolynomial) * pow( 2, tElementLevel ))
            {
                tNeighborPosition( 0 ) = tElementPosition( 0 ) + i;
                tNeighborPosition( 1 ) = tElementPosition( 1 ) + j;
                tNeighborPosition( 2 ) = tElementPosition( 2 ) + k;
                tElement = mBaseElement.give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
                if ( ElementActive.test(tElement) == 1)
                {
                    tElementNeighbor(11,tVar) = tElement;
                }
            }
            tVar++;
        }
    }
    return tElementNeighbor;
}
