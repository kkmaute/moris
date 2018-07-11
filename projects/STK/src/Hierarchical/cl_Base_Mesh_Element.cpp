/*
 * cl_Base_Mesh_Element.cpp
 *
 *  Created on: Feb 22, 2018
 *      Author: gleim
 */

#include "cl_Base_Mesh_Element.hpp"
using namespace moris;

uint
Base_Mesh_Element::give_element_of_position(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aIJKPosition) const
{
    uint tElementId = 0;
    uint tLevel = aLevel - 1;
    uint tY_Position; // Position of the Element in y-direction
    if ( aModelDim == 1 )
    {
        if ( aLevel > 0 )
        {
            //Give number of elements for the last level to have the initial number at position 0,0
            tElementId = give_number_of_elements( tLevel, aModelDim, aNumberOfElementsPerDirection );
        }
        // Count the number of elements of the lower levels , the x-position and the y-positon to get the element ID
        tElementId = tElementId + aIJKPosition( 0 );
    }
    else if ( aModelDim == 2 )
    {
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 );
        if ( aLevel > 0 )
        {
            //Give number of elements for the last level to have the initial number at position 0,0
            tElementId = give_number_of_elements( tLevel, aModelDim, aNumberOfElementsPerDirection );
            // Position of the Element in y-direction for a specific level
            tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel );
        }
        // Count the number of elements of the lower levels , the x-position and the y-positon to get the element ID
        tElementId = tElementId + aIJKPosition( 0 ) + tY_Position;
    }
    else if ( aModelDim == 3 )
    {
        //Position of the element in y- and z-direction
        tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 )
                    + aIJKPosition( 2 ) * aNumberOfElementsPerDirection( 0 ) * aNumberOfElementsPerDirection( 1 );
        if ( aLevel>0 )
        {
            //Give number of elements for the last level to have the initial number at position 0,0,0
            tElementId = give_number_of_elements( tLevel, aModelDim, aNumberOfElementsPerDirection );
            // Position of the Element in y and z--direction for a specific level
            tY_Position = aIJKPosition( 1 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel )
                        + aIJKPosition( 2 ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, aLevel )
                        * aNumberOfElementsPerDirection( 1 ) * pow( 2, aLevel );
        }
        // Count the number of elements of the lower levels , the x-position and the y- and z-positon to get the element ID
        tElementId = tElementId + aIJKPosition( 0 ) + tY_Position;
    }
    return tElementId;
}

//--------------------------------------------------------------------------------

/* @TODO give_position_of_element: cleanup formulas */
Mat<uint>
Base_Mesh_Element::give_position_of_element(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Element positions in the direction x,y,z
    Mat<uint> tElementPosition( aModelDim, 1 );
    //On which level lives this element
    uint tElementLevel=give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // temporary variable for the number of elements per level
    uint tNumberOfElements = 0;
    if ( aModelDim == 1 )
    {
        if ( tElementLevel != 0 )
        {
            // Count the elements until the level tElement_level-1
            tNumberOfElements = give_number_of_elements( tElementLevel-1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Position is id based ( starting with 1,2,...)
        tElementPosition( 0 ) = aElementId +1 - tNumberOfElements
                + pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 );
        // Position is index based ( starting with 0,1,2,...)
        tElementPosition( 0 ) -= 1;
    }
    else if ( aModelDim == 2 )
    {
        if ( tElementLevel != 0 )
        {
            // Count the elements until the level tElement_level-1
            tNumberOfElements = give_number_of_elements( tElementLevel-1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Position is id based ( starting with 1,2,...)
        tElementPosition( 1 ) = ceil( ( aElementId + 1 - tNumberOfElements )
                / ( pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 ) ) );
        // Position is id based ( starting with 1,2,...)
        tElementPosition( 0 ) = aElementId + 1 - tNumberOfElements
                + pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 )
                - tElementPosition( 1 ) * pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 );
        // Position is index based ( starting with 0,1,2,...)
        tElementPosition( 1 ) -= 1;
        tElementPosition( 0 ) -= 1;
    }
    else if ( aModelDim == 3 )
    {
        if ( tElementLevel != 0 )
        {
            // Count the elements until the level tElement_level-1
            tNumberOfElements = give_number_of_elements( tElementLevel-1, aModelDim, aNumberOfElementsPerDirection );
        }
        // Position is id based ( starting with 1,2,...)
        tElementPosition( 2 ) = ceil( ( aElementId + 1 - tNumberOfElements ) / ( pow( 2, tElementLevel )
                * aNumberOfElementsPerDirection( 0 ) * pow(2, tElementLevel) * aNumberOfElementsPerDirection( 1 )) );
        tElementPosition( 1 ) = ceil( ( aElementId + 1 - tNumberOfElements - ( tElementPosition( 2 ) - 1 ) * pow( 2, tElementLevel )
                * aNumberOfElementsPerDirection( 0 ) * pow(2, tElementLevel) * aNumberOfElementsPerDirection( 1 ) )
                / ( pow( 2, tElementLevel) * aNumberOfElementsPerDirection( 0 ) ) );
        tElementPosition( 0 ) = aElementId + 1 - tNumberOfElements + pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 )
                            - tElementPosition( 1 ) * pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 )
                            + pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 1 )
                            - tElementPosition( 2 ) * pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 0 )
                            * pow( 2, tElementLevel ) * aNumberOfElementsPerDirection( 1 );
        // Position is index based ( starting with 0,1,2,...)
        tElementPosition( 2 ) -= 1;
        tElementPosition( 1 ) -= 1;
        tElementPosition( 0 ) -= 1;
    }
    return tElementPosition;
}

//--------------------------------------------------------------------------------

/* @TODO give_element_level: cleanup formulas */
uint
Base_Mesh_Element::give_element_level(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // output variable
    uint tElementLevel = 1;
    // temporary variable for the while loop
    uint tLevel = 1;
    //Compute the relation of the different levels by the power of the level
    uint tPowLevel = pow( 2, tElementLevel - 1 );
    if ( aModelDim == 1 )
    {
        // temporary variable for the number of elements per level
        uint tNumberOfElements = tPowLevel * aNumberOfElementsPerDirection( 0 );
        while( tLevel > 0 )
        {
            // Checks, if the number of elements fits in the respective level
            tLevel = floor( (real)( aElementId + 1) / tNumberOfElements );
            if ( (real)( aElementId + 1 ) / tNumberOfElements == 1 || tLevel < 1 )
            {
                // If yes, turn the switch to stop the while loop
                tLevel = 0;
            }
            else
            {
                // If it is not in the level, add one level and add the number of elements from this level
                tElementLevel++;
                tPowLevel = pow( 2, tElementLevel - 1 );
                tNumberOfElements += tPowLevel * aNumberOfElementsPerDirection( 0 );
            }
        }
        // While loop stops always one level to late, so minus one
        tElementLevel--;
    }
    else if ( aModelDim == 2 )
    {
        // temporary variable for the number of elements per level
        uint tNumberOfElements = tPowLevel * aNumberOfElementsPerDirection( 0 )
                               * tPowLevel * aNumberOfElementsPerDirection( 1 );
        while(tLevel>0)
        {
            // Checks, if the number of elements fits in the respective level
            tLevel = floor( (real)( aElementId + 1 )/ tNumberOfElements );
            if ( (real)( aElementId + 1 )/ tNumberOfElements == 1 || tLevel < 1 )
            {
                // If yes, turn the switch to stop the while loop
                tLevel = 0;
            }
            else
            {
                // If it is not in the level, add one level and add the number of elements from this level
                tElementLevel++;
                tPowLevel = pow( 2, tElementLevel - 1 );
                tNumberOfElements += tPowLevel * aNumberOfElementsPerDirection( 0 )
                                   * tPowLevel * aNumberOfElementsPerDirection( 1 );
            }
        }
        // While loop stops always one level to late, so minus one
        tElementLevel--;
    }
    else if ( aModelDim == 3 )
    {
        // temporary variable for the number of elements per level
        uint tNumberOfElements = tPowLevel * aNumberOfElementsPerDirection( 0 )
                              * tPowLevel * aNumberOfElementsPerDirection( 1 ) * tPowLevel * aNumberOfElementsPerDirection( 2 );
        while( tLevel > 0 )
        {
            // Checks, if the number of elements fits in the respective level
            tLevel = floor( (real)( aElementId + 1 )/ tNumberOfElements );
            if ( (real)( aElementId + 1 )/ tNumberOfElements == 1 || tLevel < 1 )
            {
                // If yes, turn the switch to stop the while loop
                tLevel = 0;
            }
            else
            {
                // If it is not in the level, add one level and add the number of elements from this level
                tElementLevel++;
                tPowLevel = pow( 2, tElementLevel - 1 );
                tNumberOfElements += tPowLevel * aNumberOfElementsPerDirection( 0 )
                                   * tPowLevel * aNumberOfElementsPerDirection( 1 ) * tPowLevel * aNumberOfElementsPerDirection( 2 );
            }
        }
        tElementLevel--;
    }
    return tElementLevel;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Element::give_parent_of_element(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Computes the position of an element
    Mat<uint> tElementPosition = this->give_position_of_element( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Determines the level of the element
    uint tElementLevel = this->give_element_level( aElementId, aModelDim, aNumberOfElementsPerDirection );
    // Predefine UINT_MAX, if no parent can be found UINT_MAX will be returned
    uint tParent = UINT_MAX;
    uint tParentLevel = 0;
    if (tElementLevel>0)
    {
        // Go one level lower to check for a parent element
        tParentLevel = tElementLevel - 1;
        // Position of the Parent element
        Mat<uint> tParentPosition( aModelDim , 1 );
        if ( aModelDim == 1 )
        {
            // Two possibilities in 1D if the two children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 )  = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
            }
        }
        else if ( aModelDim == 2 )
        {
            // Four possibilities in 2D if the four children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
            }
        }
        else if ( aModelDim == 3 )
        {
            // Eight possibilities in 3D if the four children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 +1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 -1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 +1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
            }
        }
        //Determine with the position of the parent the element Id
        tParent = this->give_element_of_position( tParentLevel, aModelDim, aNumberOfElementsPerDirection, tParentPosition );
    }
    return tParent;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Element::give_parent_of_level_x(
        uint const & aElement,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aWhichLevel)
{
    // Calculates first the parent of the element
    uint tParent = this->give_parent_of_element( aElement, aModelDim, aNumberOfElementsPerDirection );
    uint tLevel = UINT_MAX;
    // If the parent is not UINT_MAX, there is a possibility to determine a parent for a specific level, otherwise return UINT_MAX
    if ( tParent < UINT_MAX)
    {
        // Parent level is no the new starting element
        tParent = aElement;
        // Determine the level for this new element
        tLevel = this->give_element_level( tParent, aModelDim, aNumberOfElementsPerDirection );
        // While loop runs, until the specific level is reached
        while( tLevel > aWhichLevel)
        {
            // Determines the parent element
            tParent = this->give_parent_of_element( tParent, aModelDim, aNumberOfElementsPerDirection );
            // Determines the element level for the while loop and the user specific level
            tLevel = this->give_element_level( tParent, aModelDim, aNumberOfElementsPerDirection );
        }
    }
    return tParent;
}

Mat<uint>
Base_Mesh_Element::give_parent_child_realation_of_element(
        uint const & aElement,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection)
{
    // Computes the position of an element
    Mat<uint> tElementPosition = this->give_position_of_element( aElement, aModelDim, aNumberOfElementsPerDirection );
    // Determines the level of the element
    uint tElementLevel = this->give_element_level( aElement, aModelDim, aNumberOfElementsPerDirection );
    // Predefine UINT_MAX, Matrix with first entry = parent and second entry the child relation
    Mat<uint> tParentChild( 2 , 1 , UINT_MAX );
    uint tParentLevel = 0;
    if ( tElementLevel > 0 )
    {
        // Go one level lower to check for a parent element
        tParentLevel = tElementLevel - 1;
        // Position of the Parent element
        Mat<uint> tParentPosition( aModelDim , 1 );
        if ( aModelDim == 1 )
        {
            // Two possibilities in 1D if the two children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 )  % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 0;
            }
        }
        else if ( aModelDim == 2 )
        {
            // Four possibilities in 2D if the four children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 3;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 2;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 0;
            }
        }
        else if ( aModelDim == 3 )
        {
            // Eight possibilities in 3D if the four children have even or odd position
            if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 7;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 5;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 6;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) == 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 ) / 2 - 1;
                tParentChild( 1 ) = 4;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 3;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 1;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) == 0
                    && ( ( tElementPosition( 2 ) +1 ) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 2;
            }
            else if ( ( ( tElementPosition( 0 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 1 ) + 1 ) % 2 ) != 0
                    && ( ( tElementPosition( 2 ) + 1) % 2 ) != 0 )
            {
                tParentPosition( 0 ) = ( tElementPosition( 0 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 1 ) = ( tElementPosition( 1 ) + 1 + 1 ) / 2 - 1;
                tParentPosition( 2 ) = ( tElementPosition( 2 ) + 1 + 1 ) / 2 - 1;
                tParentChild( 1 ) = 0;
            }
        }
        //Determine with the position of the parent the element Id
        tParentChild( 0 ) = this->give_element_of_position( tParentLevel, aModelDim, aNumberOfElementsPerDirection, tParentPosition );
    }
    return tParentChild;
}

//--------------------------------------------------------------------------------
/* @TODO give_number_of_elements: cleanup formulas */
uint
Base_Mesh_Element::give_number_of_elements(
        uint const & aLevel,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection)
{
    // output variable
    uint tElementNumber = 0;
    if ( aModelDim == 1 )
    {
        for ( uint i = 0; i < aLevel + 1; i++ )
        {
            // Loop over the user specific level and count all elements from each level
            tElementNumber += pow( 2, i ) * aNumberOfElementsPerDirection( 0 );
        }
    }
    else if ( aModelDim == 2 )
    {
        for ( uint i = 0; i < aLevel + 1; i++ )
        {
            // Loop over the user specific level and count all elements from each level
            tElementNumber += pow( 2, i ) * aNumberOfElementsPerDirection( 0 )
                           * pow( 2, i ) * aNumberOfElementsPerDirection( 1 );
        }
    }
    else if ( aModelDim == 3 )
    {
        for ( uint i = 0; i < aLevel + 1; i++ )
        {
            // Loop over the user specific level and count all elements from each level
            tElementNumber += pow( 2, i ) * aNumberOfElementsPerDirection( 0 )
                           * pow( 2, i ) * aNumberOfElementsPerDirection( 1 )
                           * pow( 2, i ) * aNumberOfElementsPerDirection( 2 );
        }
    }
    return tElementNumber;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Element::give_children_of_element(
        uint const & aElement,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Determines the position of the element
    Mat<uint> tElementPosition = this->give_position_of_element( aElement, aModelDim, aNumberOfElementsPerDirection );
    // Determines the level of the element and add one to get the level of the children level
    uint tChildrenLevel = this->give_element_level( aElement, aModelDim, aNumberOfElementsPerDirection ) + 1;
    // Vector with all children of the element aElement
    Mat<uint> tChildren( pow( 2, aModelDim ) , 1 );
    // Position of the children with respect to the element aElement
    Mat<uint> tChildrenPosition( aModelDim , 1 );
    if ( aModelDim == 1 )
    {
        // Determines the position of the children
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        // Computes the element Id from the position
        tChildren( 0 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildren( 1 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
    }
    else if ( aModelDim == 2 )
    {
        // Determines the position of the children
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        // Computes the element Id from the position
        tChildren( 0 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        tChildren( 1 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildren( 2 ) =this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildren( 3 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
    }
    else if ( aModelDim == 3 )
    {
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1 - 1;
        tChildren( 0 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1 - 1;
        tChildren( 1 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1 - 1;
        tChildren( 2 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1 - 1;
        tChildren( 3 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1;
        tChildren( 4 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1 - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1;
        tChildren( 5 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1 - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1;
        tChildren( 6 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
        tChildrenPosition( 0 ) = 2 * ( tElementPosition( 0 ) + 1 ) - 1;
        tChildrenPosition( 1 ) = 2 * ( tElementPosition( 1 ) + 1 ) - 1;
        tChildrenPosition( 2 ) = 2 * ( tElementPosition( 2 ) + 1 ) - 1;
        tChildren( 7 ) = this->give_element_of_position( tChildrenLevel, aModelDim, aNumberOfElementsPerDirection, tChildrenPosition );
    }
    return tChildren;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Element::give_neighbor_of_element(
        uint const & aElement,
        uint const & aModelDim,
        uint const & aBuffer,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    // Deterimines the position of the element
    Mat<uint> tElementPosition = give_position_of_element( aElement, aModelDim, aNumberOfElementsPerDirection );
    // Computes the level of the element
    uint tElementLevel = give_element_level( aElement, aModelDim, aNumberOfElementsPerDirection );
    // Element plus neighbors (9 for aModelDim = 2, 27 for aModelDim = 3) if a buffer of one, otherwise more
    Mat<uint> tElementNeighbor( 1 , pow( 1 + 2 * aBuffer, aModelDim ) , UINT_MAX );
    // Position of a neighbor element;
    Mat<uint> tNeighborPosition( aModelDim , 1 );
    if ( aModelDim == 1 )
    {
        if ( aBuffer == 1 )
        {
            // Determine the position of the neighbor with respect to the position of the element
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            // Determine the element id of the neighbor
            tElementNeighbor( 0 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tElementNeighbor( 1 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tElementNeighbor( 2 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else if ( aBuffer > 1)
        {
            // If a higher buffer layer is needed, a loop generates the neighbors
            uint tVar = 0;
            Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
            for ( int i = -( (int)aBuffer ); i <= ( (int)aBuffer ); i++ )
            {
                if ( tElementPosition( 0 ) + i >= 1 * pow( 2, tElementLevel )
                && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel ) )
                {
                    aIJKPosition( 0 ) = tElementPosition( 0 ) + i;
                    tElementNeighbor(tVar) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, aIJKPosition );
                }
                tVar++;
            }
        }
    }
    else if ( aModelDim == 2 )
    {
        if ( aBuffer == 1 )
        {
            // Determine the position of the neighbor with respect to the position of the element
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            // Determine the element id of the neighbor
            tElementNeighbor( 0 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tElementNeighbor( 1 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tElementNeighbor( 2 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tElementNeighbor( 3 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tElementNeighbor( 4 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tElementNeighbor( 5 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tElementNeighbor( 6 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tElementNeighbor( 7 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tElementNeighbor( 8 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else if ( aBuffer > 1)
        {
            // If a higher buffer layer is needed, a loop generates the neighbors
            uint tVar = 0;
            Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
            for ( int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
            {
                for ( int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
                {
                    if ( tElementPosition( 0 ) + i >= 1 * pow( 2, tElementLevel )
                    && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel )
                    && tElementPosition( 1 ) + j >= 1 * pow( 2, tElementLevel )
                    && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel ) )
                    {
                        aIJKPosition( 0 ) = tElementPosition( 0 ) + i;    aIJKPosition( 1 ) = tElementPosition( 1 ) + j;
                        tElementNeighbor(tVar) = this->give_element_of_position(tElementLevel,aModelDim,aNumberOfElementsPerDirection,aIJKPosition);
                    }
                    tVar++;
                }
            }
        }
    }
    else if (aModelDim == 3)
    {
        if ( aBuffer == 1)
        {
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 0 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 1 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 2 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 3 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 4 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 5 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 6 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 7 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) - 1;
            tElementNeighbor( 8 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 9 )  = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 10 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 11 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 12 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 13 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 14 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 15 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 16 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 );
            tElementNeighbor( 17 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 18 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 19 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) - 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 20 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 21 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 22 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 );
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 23 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) - 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 24 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 );
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor( 25 ) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
            tNeighborPosition( 0 ) = tElementPosition( 0 ) + 1;
            tNeighborPosition( 1 ) = tElementPosition( 1 ) + 1;
            tNeighborPosition( 2 ) = tElementPosition( 2 ) + 1;
            tElementNeighbor(26) = this->give_element_of_position( tElementLevel, aModelDim, aNumberOfElementsPerDirection, tNeighborPosition );
        }
        else if ( aBuffer > 1)
        {
            // If a higher buffer layer is needed, a loop generates the neighbors
            uint tVar = 0;
            Mat<uint> aIJKPosition( aModelDim , 1 , 0 );
            for ( int k = -((int)aBuffer); k <= ((int)aBuffer); k++ )
            {
                for ( int j = -((int)aBuffer); j <= ((int)aBuffer); j++ )
                {
                    for ( int i = -((int)aBuffer); i <= ((int)aBuffer); i++ )
                    {
                        if ( tElementPosition( 0 ) + i >= 1 * pow( 2, tElementLevel )
                        && tElementPosition( 0 ) + i < aNumberOfElementsPerDirection( 0 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel ) &&
                        tElementPosition( 1 ) + j >= 1 * pow( 2, tElementLevel )
                        && tElementPosition( 1 ) + j < aNumberOfElementsPerDirection( 1 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel ) &&
                        tElementPosition( 2 ) + k >= 1 * pow( 2, tElementLevel )
                        && tElementPosition( 2 ) + k < aNumberOfElementsPerDirection( 2 ) * pow( 2, tElementLevel ) - 1 * pow( 2, tElementLevel ) )
                        {
                            aIJKPosition( 0 ) = tElementPosition( 0 ) + i;
                            aIJKPosition( 1 ) = tElementPosition( 1 ) + j;
                            aIJKPosition( 2 ) = tElementPosition( 2 ) + k;
                            tElementNeighbor(tVar) = this->give_element_of_position(tElementLevel,aModelDim,aNumberOfElementsPerDirection,aIJKPosition);
                        }
                        tVar++;
                    }
                }
            }
        }
    }
    return tElementNeighbor;
}

//--------------------------------------------------------------------------------

uint
Base_Mesh_Element::give_element_owner(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> & aProcRange,
        Mat<uint> const & aProcneighbors)
{
    // Determines the position of the element
    Mat<uint> tElementPosition = this->give_position_of_element(aElementId,aModelDim,aNumberOfElementsPerDirection);
    // Computes the level of the element
    uint tElementLevel = this->give_element_level(aElementId,aModelDim,aNumberOfElementsPerDirection);
    uint tProcOwner = UINT_MAX;
    if ( par_size() > 1 )
    {
        uint tProcRank = par_rank();
        // Calculate the range of the proc for the level of the element
        Mat<uint> tProcRange( aProcRange.length() , 1 , 0 );
        // X-Position start of the proc
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tElementLevel );
        // X-Position end of the proc
        tProcRange( 1 ) = ( aProcRange( 1 ) + 1 ) * pow( 2, tElementLevel ) - 1;
        // Y-Position start of the proc
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tElementLevel );
        // Y-Position end of the proc
        tProcRange( 3 ) = ( aProcRange( 3 ) + 1 ) * pow( 2, tElementLevel ) - 1;
        //Different possibilities for an ownership, starting from the highest possible neighbor (see Hierarchical_Mesh_Main::proc_coordinates_neighbors())
        MORIS_ASSERT( aModelDim>1,"Element ownership is not implemented for 1D");
        if ( aModelDim == 2 )
        {
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 ) )
                tProcOwner =tProcRank;
            else if ( tElementPosition( 1 ) == tProcRange( 2 ) - 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 ) )
                tProcOwner = aProcneighbors( 0 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 ) )
                tProcOwner = aProcneighbors( 1 );
            else if ( tElementPosition( 1 ) == tProcRange( 3 ) + 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 ) )
                tProcOwner = aProcneighbors( 2 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 ) )
                tProcOwner = aProcneighbors( 3 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1 )
                tProcOwner = aProcneighbors( 6 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1 )
                tProcOwner = aProcneighbors( 15 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1 )
                tProcOwner = aProcneighbors( 16 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1 )
                tProcOwner = aProcneighbors( 17 );
        }
        else if ( aModelDim == 3)
        {
            // Z-Position start of the proc
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tElementLevel );
            // Z-Position end of the proc
            tProcRange( 5 ) = ( aProcRange( 5 ) + 1 ) * pow( 2, tElementLevel ) - 1;
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner =tProcRank;
            else if ( tElementPosition( 1 ) == tProcRange( 2 ) - 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 ) && tElementPosition( 2 ) >= tProcRange( 4 )
                    && tElementPosition( 2 ) <= tProcRange( 5 )  )
                tProcOwner = aProcneighbors( 0 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 )  && tElementPosition( 2 ) >= tProcRange( 4 )
                    && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 1 );
            else if ( tElementPosition( 1 ) == tProcRange( 3 ) + 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 )  && tElementPosition( 2 ) >= tProcRange( 4 )
                    && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 2 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 )  && tElementPosition( 2 ) >= tProcRange( 4 )
                    && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 3 );
            else if ( tElementPosition( 2 ) == tProcRange( 4 ) - 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 ) && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 )  )
                tProcOwner = aProcneighbors( 4 );
            else if ( tElementPosition( 2 ) == tProcRange( 5 ) + 1 && tElementPosition( 0 ) >= tProcRange( 0 )
                    && tElementPosition( 0 ) <= tProcRange( 1 ) && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 )  )
                tProcOwner = aProcneighbors( 5 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 8 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 15 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 16 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
                tProcOwner = aProcneighbors( 17 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 6 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 10 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 12 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 14 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 18 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 2 ) - 1
                    && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 20 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 23 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) == tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 25 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 ) && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 7 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 ) && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 11 );
            else if ( tElementPosition( 0 ) == tProcRange( 0 ) - 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 ) && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 21 );
            else if ( tElementPosition( 0 ) == tProcRange( 1 ) + 1 && tElementPosition( 1 ) >= tProcRange( 2 )
                    && tElementPosition( 1 ) <= tProcRange( 3 )&& tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 22 );
            else if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) == tProcRange( 2 ) - 1 && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 9 );
            else if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) == tProcRange( 3 ) + 1 && tElementPosition( 2 ) == tProcRange( 4 ) - 1  )
                tProcOwner = aProcneighbors( 13 );
            else if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) == tProcRange( 2 ) - 1 && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 19 );
            else if ( tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) == tProcRange( 3 ) + 1 && tElementPosition( 2 ) == tProcRange( 5 ) + 1  )
                tProcOwner = aProcneighbors( 24 );
        }
    }
    return tProcOwner;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Element::give_element_share(
        uint const & aElementId,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> & aProcRange,
        Mat<uint> const & aProcNeighbors)
{
    // Computes the position of the element
    Mat<uint> tElementPosition = this->give_position_of_element(aElementId,aModelDim,aNumberOfElementsPerDirection);
    // Determines the level of the element
    uint tElementLevel = this->give_element_level(aElementId,aModelDim,aNumberOfElementsPerDirection);
    // Arbitrary number of possible neighbors, Proc numbers will be found more then once, because of different criterias, Unique defines then a unique list of processors
    Mat<uint> tProcShare( 40 , 1 , UINT_MAX );
    uint tVar = 0;
    MORIS_ASSERT( aModelDim>1,"Element share is not implemented for 1D");
    if ( par_size() > 1 )
    {
        uint tProcRank = par_rank();
        // Calculate the range of the proc for the level of the element
        Mat<uint> tProcRange( aProcRange.length() , 1 , 0 );
        // X-Position start of the proc
        tProcRange( 0 ) = aProcRange( 0 ) * pow( 2, tElementLevel );
        // X-Position end of the proc
        tProcRange( 1 ) = ( aProcRange( 1 ) + 1 ) * pow( 2, tElementLevel ) - 1;
        // Y-Position start of the proc
        tProcRange( 2 ) = aProcRange( 2 ) * pow( 2, tElementLevel );
        // Y-Position end of the proc
        tProcRange( 3 ) = ( aProcRange( 3 ) + 1 ) * pow( 2, tElementLevel ) - 1;
        //Checks for different possible scenarios (Element is along an edge or in a corner of a proc or in the aura of a proc
        if ( aModelDim == 2 )
        {
            if ( tElementPosition( 0 ) >= (tProcRange( 0 ) - 1) && tElementPosition( 0 ) <= (tProcRange( 1 ) + 1)
                    && tElementPosition( 1 ) >= (tProcRange( 2 ) - 1) && tElementPosition( 1 ) <= (tProcRange( 3 ) + 1) )
            {
                tProcShare(tVar) = tProcRank;
                tVar++;
            }
            if ( tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= (tProcRange( 1 )) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
            }
            if ( tElementPosition( 0 ) <= tProcRange( 1 ) + 1 && tElementPosition( 0 ) >= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= (tProcRange( 3 ))  )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
            }
            if ( tElementPosition( 1 ) <= tProcRange( 3 ) + 1 && tElementPosition( 1 ) >= tProcRange( 3 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= (tProcRange( 1 ) ) )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= (tProcRange( 3 ) )  )
            {
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= (tProcRange( 1 ) + 1)
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 15 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= (tProcRange( 1 ) + 1)
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= tProcRange( 3 ) + 1   )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 17 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= (tProcRange( 0 ))
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= tProcRange( 3 ) + 1   )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 16 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= (tProcRange( 0 ))
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )   )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 6 );
                tVar++;
            }
        }
        else if ( aModelDim == 3 )
        {
            // Z-Position start of the proc
            tProcRange( 4 ) = aProcRange( 4 ) * pow( 2, tElementLevel );
            // Z-Position end of the proc
            tProcRange( 5 ) = ( aProcRange( 5 ) + 1 ) * pow( 2, tElementLevel ) - 1;
            if ( tElementPosition( 0 ) >= (tProcRange( 0 ) - 1) && tElementPosition( 0 ) <= (tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= (tProcRange( 2 ) - 1) && tElementPosition( 1 ) <= ( tProcRange( 3 ) + 1 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 ) - 1) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = tProcRank;
                tVar++;
            }
            if ( tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 )) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
            }
            if ( tElementPosition( 0 ) <= tProcRange( 1 ) + 1 && tElementPosition( 0 ) >= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= (tProcRange( 3 ))
                    && tElementPosition( 2 ) >= (tProcRange( 4 )) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
            }
            if ( tElementPosition( 1 ) <= tProcRange( 3 ) + 1 && tElementPosition( 1 ) >= tProcRange( 3 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ))
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= ( tProcRange( 1 ) + 1)
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 )) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 15 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= (tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 17 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= tProcRange( 3 ) + 1
                    && tElementPosition( 2 ) >= (tProcRange( 4 )) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 16 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= tProcRange( 4 ) && tElementPosition( 2 ) <= tProcRange( 5 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 8 );
                tVar++;
            }
            if ( tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 ) - 1) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 9 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) <= tProcRange( 1 ) + 1 && tElementPosition( 0 ) >= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 11 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 1 ) <= tProcRange( 3 ) + 1 && tElementPosition( 1 ) >= tProcRange( 3 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 13 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 7 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= ( tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 ) - 1) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 15 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 9 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 10 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 11 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= ( tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= ( tProcRange( 3 ) + 1 )
                    && tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 17 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 11 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 13 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 14 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= ( tProcRange( 3 ) + 1 )
                    && tElementPosition( 2 ) >= (tProcRange( 4 ) - 1) && tElementPosition( 2 ) <= tProcRange( 4 )  )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 16 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 13 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 7 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 12 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 )  )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 8 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 9 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 7 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 6 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
            if ( tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= (tProcRange( 1 ))
                    && tElementPosition( 2 ) >= (tProcRange( 5 )) && tElementPosition( 2 ) <= (tProcRange( 5 ) + 1) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 19 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) <= tProcRange( 1 ) + 1 && tElementPosition( 0 ) >= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 22 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 1 ) <= tProcRange( 3 ) + 1 && tElementPosition( 1 ) >= tProcRange( 3 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 24 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 21 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= ( tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 15 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 19 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 20 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 22 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 1 ) && tElementPosition( 0 ) <= ( tProcRange( 1 ) + 1 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= ( tProcRange( 3 ) + 1 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 1 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 17 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 22 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 24 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 25 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 3 ) && tElementPosition( 1 ) <= ( tProcRange( 3 ) + 1 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 2 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 16 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 21 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 23 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 24 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 0 ) >= tProcRange( 0 ) - 1 && tElementPosition( 0 ) <= tProcRange( 0 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) - 1 && tElementPosition( 1 ) <= tProcRange( 2 )
                    && tElementPosition( 2 ) >= tProcRange( 5 ) && tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 0 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 3 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 8 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 18 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 19 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 21 );
                tVar++;
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 2 ) <= ( tProcRange( 5 ) + 1 ) && tElementPosition( 2 ) >= tProcRange( 5 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 5 );
                tVar++;
            }
            if ( tElementPosition( 2 ) >= ( tProcRange( 4 ) - 1 ) && tElementPosition( 2 ) <= tProcRange( 4 )
                    && tElementPosition( 0 ) >= tProcRange( 0 ) && tElementPosition( 0 ) <= tProcRange( 1 )
                    && tElementPosition( 1 ) >= tProcRange( 2 ) && tElementPosition( 1 ) <= tProcRange( 3 ) )
            {
                tProcShare(tVar) = aProcNeighbors( 4 );
                tVar++;
            }
        }
    }
    //Resize to the found number of processors, which share this element
    tProcShare.resize( tVar, 1 );
    // Take only a unique list of processors, mulitple numbers are not needed
    tProcShare = unique( tProcShare );
    // If the proc is not owner or not sharing or doesnt no it (these two criterias), then resize to zero
    if ( isempty( tProcShare ) == 0 && tProcShare(tProcShare.length() - 1 ) == UINT_MAX )
    {
        tProcShare.resize( tProcShare.length() - 1, 1 );
    }
    return tProcShare;
}

//--------------------------------------------------------------------------------

Mat<uint>
Base_Mesh_Element::give_face_neighbor(
        uint const & aElement,
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection) const
{
    //Neighbors of the element
    Mat<uint> tNeighborElements = this->give_neighbor_of_element( aElement, aModelDim, 1, aNumberOfElementsPerDirection );
    // Only direct neigbors, which are connected to an edge or face, no diagonal neighbors
    Mat<uint> tNeighborElementsMod( 2 * aModelDim , 1 , UINT_MAX );
    MORIS_ASSERT( aModelDim>1,"Floodfill is not implemented for 1D");
    if ( aModelDim == 2 )
    {
        tNeighborElementsMod( 0 ) = tNeighborElements( 1 ); //Bottom neighbor
        tNeighborElementsMod( 1 ) = tNeighborElements( 3 ); //Left neighbor
        tNeighborElementsMod( 2 ) = tNeighborElements( 5 ); //Right neighbor
        tNeighborElementsMod( 3 ) = tNeighborElements( 7 ); //Top neighbor
    }
    else if ( aModelDim == 3)
    {
        tNeighborElementsMod( 0 ) = tNeighborElements( 10 ); //Bottom neighbor
        tNeighborElementsMod( 1 ) = tNeighborElements( 12 ); //Left neighbor
        tNeighborElementsMod( 2 ) = tNeighborElements( 14 ); //Right neighbor
        tNeighborElementsMod( 3 ) = tNeighborElements( 16 ); //Top neighbor
        tNeighborElementsMod( 4 ) = tNeighborElements( 4 );  //Front neighbor
        tNeighborElementsMod( 5 ) = tNeighborElements( 22 ); //Back neighbor
    }
    return tNeighborElementsMod;
}
