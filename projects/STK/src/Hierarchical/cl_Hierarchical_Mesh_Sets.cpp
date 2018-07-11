/*
 * cl_Hierarchical_Mesh_Sets.cpp
 *
 *  Created on: Jan 5, 2018
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Sets.hpp" // STK/src/Hierarchical
using namespace moris;
Mat<uint>
Hierarchical_Mesh_Sets::set_nodeset(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aInNodeSet,
        Mat<uint> const & aElementListOnProc)
{
    //Create a node set
    Mat<uint> tElementOfBasis;
    uint tVar = 0;
    //Three Numbers per node are needed (For Example Level  = 0, Node = 123, Node set = 1 )
    Mat<uint> aOutNodeSet( aInNodeSet.n_rows(), 3, 0 );
    Mat<uint> tHelpMat;
    uint tHelp;
    for ( uint i=0; i < aInNodeSet.n_rows(); i++ )
    {
        tElementOfBasis = mBasis.give_element_of_basis( aInNodeSet( i, 1 ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        for( uint j = 0; j < tElementOfBasis.length(); j++ )
        {
            tHelpMat = ( tElementOfBasis(j) == aElementListOnProc );
            tHelp = sum( tHelpMat );
            if ( tHelp > 0 )
            {
                // Level  = 0, Node = 1, Node set = 2
                aOutNodeSet( tVar, 0 ) = aInNodeSet( i, 0 );
                aOutNodeSet( tVar, 1 ) = aInNodeSet( i, 1 );
                aOutNodeSet( tVar, 2 ) = aInNodeSet( i, 2 );
                tVar++;
                break;
            }
        }
    }
    aOutNodeSet.resize( tVar, 3 );
    return aOutNodeSet;
}

Cell<Mat<uint>>
Hierarchical_Mesh_Sets::update_nodeset(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aLevel,
        Mat<uint> const & aNumSets,
        Mat<uint> & aNodeSet,
        BoostBitset const & aBasisActive)
{
    // Initilize Cell for mesh output
    Cell<Mat<uint>> tSetsEntIds( aNumSets( 0 ) );
    uint tVarCell = 0;
    MORIS_ASSERT( aPolynomial == 1, "Only for a linear polynomial possible, otherwise the B-spline coefficients do not overlapp!" );
    if ( aPolynomial == 1 )
    {
        if ( aNumSets( 0 ) > 0 )
        {
            if ( aNodeSet.n_rows() > 0)
            {
                // Update NodeSet
                // Which max level is in the NodeSet list
                uint tHelp = aNodeSet( aNodeSet.n_rows() - 1, 0 );
                // How many NodeSets exsists
                uint Length = aNodeSet.n_rows();
                if ( tHelp < aLevel )
                {
                    aNodeSet.resize( aNodeSet.n_rows() * ( aLevel + 1 ), 3 );
                    Mat<uint> tElement_of_basis( pow( aPolynomial + 1, aModelDim ), 1 );
                    Mat<uint> tChildren( pow( 2, aModelDim ), 1 );
                    Mat<uint> tBasis_of_children( pow( aPolynomial + 1, aModelDim ), 1 );
                    Mat<uint> tBasis_of_element( pow( aPolynomial + 1, aModelDim ), 1 );
                    // Temporary variable for loop
                    Mat<uint> tHelpMat;
                    // Temporary variable for loop
                    uint WhichBasis;
                    // Temporary variable for loop
                    uint tVar = Length;
                    for ( uint i = 0; i < Length * aLevel; i++ )
                    {
                        WhichBasis = aNodeSet( i, 1 );
                        //Compute the elements, in which the basis has support
                        tElement_of_basis = mBasis.give_element_of_basis( WhichBasis, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                        //Compute children of first children
                        tChildren = mBaseElement.give_children_of_element( tElement_of_basis(0), aModelDim, aNumberOfElementsPerDirection );
                        //Compute basis functions of children element
                        tBasis_of_element= mHMRElement.give_basis_of_element( tElement_of_basis(0), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                        //Compare basis functions
                        tHelpMat = ( tBasis_of_element == WhichBasis );
                        tHelpMat = find( tHelpMat );
                        tHelp = tChildren( tHelpMat( 0 ) );
                        tBasis_of_children = mHMRElement.give_basis_of_element( tHelp, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
                        tHelp = tBasis_of_children( tHelpMat( 0 ) );
                        aNodeSet( tVar, 0 ) = aNodeSet( i, 0 ) + 1;
                        aNodeSet( tVar, 1 ) = tHelp;
                        aNodeSet( tVar, 2 ) = aNodeSet( i, 2 );
                        tVar++;
                    }
                    aNodeSet.resize( tVar, 3 );
                }
            }
            Mat<uint> tMatNodeSet;
            for ( uint i = 0; i < aNumSets( 0 ); i++ )
            {
                uint tVar = 0;
                if ( aNodeSet.n_rows() > 0 )
                {
                    tMatNodeSet.set_size( aNodeSet.n_rows(), 1 );
                    for( uint j = 0; j < aNodeSet.n_rows(); j++ )
                    {
                        if ( aBasisActive.test( aNodeSet( j, 1 ) ) == 1 && aNodeSet( j, 2 ) == i + 1 )
                        {
                            tMatNodeSet( tVar ) = aNodeSet( j, 1 ) +  1;
                            tVar++;
                        }
                    }
                }
                if ( tVar > 0 )
                {
                    tMatNodeSet.resize( tVar, 1 );
                    tSetsEntIds( tVarCell ) = tMatNodeSet;
                    tVarCell++;
                }
                else
                {
                    tVarCell++;
                }
            }
        }
    }
    return tSetsEntIds;
}

Mat<uint>
Hierarchical_Mesh_Sets::set_sideset(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aDecomp,
        Mat<uint> const & aSideSet)
{
    // Side set for the elements
    Mat<uint> tSideSet( ( aDecomp( 1 ) - aDecomp( 0 ) + 1 ) * ( aDecomp( 3 ) - aDecomp( 2 ) + 1 ) * 2
            + ( aDecomp( 1 ) - aDecomp( 0 ) + 1 ) * ( aDecomp( 5 ) - aDecomp( 4 ) + 1 ) * 2
            + ( aDecomp( 3 ) - aDecomp( 2 ) + 1 ) * ( aDecomp( 5 ) - aDecomp( 4 ) + 1 ) * 2, 3, 0 );
    // Element:
    //             SideSet(3)
    //            -------------
    //            |           |
    // SideSet(4) |  Element  | SideSet(2)
    //            |           |
    //            -------------
    //              SideSet(1)

    //Elements in the second and second last row get a side set on the bottom and top, respectively.
    //Temporary variable for a loop
    uint j;
    //Temporary variable for a loop
    uint count = 0;
    // Determine the initial mesh for level zero
    uint tLevel = 0;
    Mat<uint> tPosition( 3, 1 );
    for ( uint k = aDecomp( 4 ); k <= aDecomp( 5 ); k++ )
    {
        for ( uint i = aDecomp( 0 ); i <= aDecomp( 1 ); i++ )
        {
            //Check Side set in y-direction
            if ( aDecomp( 2 ) == aPolynomial)
            {
                // Starting in the second row, since the first row is a customized aura
                j=aPolynomial;
                // tPosition of the element in a tensorial grid
                tPosition( 0 ) = i;
                tPosition( 1 ) = j;
                tPosition( 2 ) = k;
                //This value means, that the side Set is at the bottom of the element
                tSideSet( count, 0 ) = 1;
                // Which element gets the side set
                tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                // A number must be chosen, for example 1 for side set 1
                tSideSet( count, 2 ) = aSideSet( 0 );
                count++;
            }
            if ( aDecomp( 3 ) == aNumberOfElementsPerDirection( 1 ) - ( aPolynomial + 1 ) )
            {
                // Starting in the second last row, since the last row is a customized aura
                j = aNumberOfElementsPerDirection( 1 )-( aPolynomial + 1 );
                // tPosition of the element in a tensorial grid
                tPosition(0) = i;
                tPosition(1) = j;
                tPosition(2) = k;
                //This value means, that the side Set is at the top of the element
                tSideSet( count, 0 ) = 3;
                // Which element gets the side set
                tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                // A number must be chosen, for example 1 for side set 1
                tSideSet( count, 2 ) = aSideSet( 2 );
                count++;
            }
        }
    }
    //Elements in the second and second last col get a side set on the left and right, respectively.
    uint i; //Temporary variable for a loop
    for ( uint k = aDecomp( 4 ); k <= aDecomp( 5 ); k++ )
    {
        for ( uint j = aDecomp( 2 ); j <= aDecomp( 3 ); j++ )
        {
            //Check Side set in x-direction
            if ( aDecomp( 0 ) == aPolynomial )
            {
                // Starting in the second col, since the first col is a customized aura
                i = aPolynomial;
                // tPosition of the element in a tensorial grid
                tPosition( 0 ) = i;
                tPosition( 1 ) = j;
                tPosition( 2 ) = k;
                //This value means, that the side Set is at the left side of the element
                tSideSet( count, 0 ) = 4;
                // Which element gets the side set
                tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                // A number must be chosen, for example 1 for side set 1
                tSideSet( count, 2 ) = aSideSet( 3 );
                count++;
            }
            if ( aDecomp( 1 ) == aNumberOfElementsPerDirection(0) - ( aPolynomial + 1 ) )
            {
                // Starting in the second last col, since the last col is a customized aura
                i = aNumberOfElementsPerDirection( 0 ) - ( aPolynomial + 1 );
                // tPosition of the element in a tensorial grid
                tPosition( 0 ) = i;
                tPosition( 1 ) = j;
                tPosition( 2 ) = k;
                //This value means, that the side Set is at the right side of the element
                tSideSet( count, 0 ) = 2;
                // Which element gets the side set
                tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                // A number must be chosen, for example 1 for side set 1
                tSideSet( count, 2 ) = aSideSet( 1 );
                count++;
            }
        }
    }
    //Elements in the second and second last col get a side set on the left and right, respectively.
    if ( aModelDim == 3 )
    {
        //Temporary variable for a loop
        uint k;
        for ( uint j = aDecomp( 2 ); j <= aDecomp( 3 ); j++ )
        {
            for ( uint i = aDecomp( 0 ); i <= aDecomp( 1 ); i++ )
            {
                //Check Side set in x-direction
                if ( aDecomp( 4 ) == aPolynomial )
                {
                    // Starting in the second col, since the first col is a customized aura
                    k = aPolynomial;
                    // tPosition of the element in a tensorial grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;
                    tPosition( 2 ) = k;
                    //This value means, that the side Set is at the left side of the element
                    tSideSet( count, 0 ) = 5;
                    // Which element gets the side set
                    tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                    // A number must be chosen, for example 1 for side set 1
                    tSideSet( count, 2 ) = aSideSet( 4 );
                    count++;
                }
                if ( aDecomp( 5 ) == aNumberOfElementsPerDirection( 2 ) - ( aPolynomial + 1 ) )
                {
                    // Starting in the second last col, since the last col is a customized aura
                    k = aNumberOfElementsPerDirection( 2 ) - ( aPolynomial + 1 );
                    // tPosition of the element in a tensorial grid
                    tPosition( 0 ) = i;
                    tPosition( 1 ) = j;
                    tPosition( 2 ) = k;
                    //This value means, that the side Set is at the right side of the element
                    tSideSet( count, 0 ) = 6;
                    // Which element gets the side set
                    tSideSet( count, 1 ) = mBaseElement.give_element_of_position( tLevel, aModelDim, aNumberOfElementsPerDirection, tPosition );
                    // A number must be chosen, for example 1 for side set 1
                    tSideSet( count, 2 ) = aSideSet( 5 );
                    count++;
                }
            }
        }
    }
    //Put the side sets in a data structure with all the importent information (see below)
    Mat<uint> tHelpMat = tSideSet.col( 2 );
    tHelpMat = tHelpMat != 0;
    uint tHelp = sum( tHelpMat );
    tHelpMat = find( tHelpMat );
    // Each element has always for side sets, which can be also zero.
    Mat<uint> tOutSideSet( tHelp, 8, 0 );
    for ( uint i = 0; i < tHelp; i++ )
    {
        // Level of the Element with side set
        tOutSideSet( i, 0 ) = 0;
        // Element number with side set
        tOutSideSet( i, 1 ) = tSideSet( tHelpMat( i ), 1 );
        // Side sets of the element
        tOutSideSet( i, 2 + tSideSet( tHelpMat( i ), 0 ) - 1 ) = tSideSet( tHelpMat( i ), 2 );
    }
    return tOutSideSet;
}

Cell<Mat<uint>>
Hierarchical_Mesh_Sets::update_sideset(
        const uint        & aModelDim,
        const Mat<uint>   & aNumberOfElementsPerDirection,
        const uint        & aLevel,
        const Mat<uint>   & aNumSets,
        const Mat<uint>   & aSideSet,
        const BoostBitset & aElementActive)
{
    // Initilize Cell for mesh output
    Cell<Mat<uint>> tSetsEntIds( aNumSets( 1 ) );
    uint tVarCell = 0;
    //Update SideSet
    Mat<uint> tSideSet = aSideSet;
    if ( aNumSets( 1 ) > 0 )
    {
        if ( aSideSet.n_rows() > 0 )
        {
            // Which max level is in the SideSet list
            uint tHelp = aSideSet( aSideSet.n_rows() - 1, 0 );
            // How many Side sets exsists
            uint Length = aSideSet.n_rows();
            // Temporary variable for loop
            uint tVar = Length;
            Mat<uint> tChildren( 1, pow( 2, aModelDim ), 0 );
            if ( tHelp < aLevel )
            {
                if ( aModelDim == 2 )
                {
                    tHelp = 0;
                    for ( uint level = 0; level <= aLevel; level++ )
                    {
                        // 1 Sideset is always divided in 2 side sets in the children
                        tHelp += pow( 2, level );
                    }
                    // Level, Element, 4 Side sets = 6 cols
                    tSideSet.resize( Length * tHelp, 6 );
                    for ( uint i = 0; i < Length * ( tHelp - pow( 2, aLevel ) ); i++ )
                    {
                        if ( tSideSet( i, 0 ) >= aLevel )
                        {
                            //If finest level of updated side sets is reached, stop that loop
                            break;
                        }
                        tChildren = mBaseElement.give_children_of_element( tSideSet( i, 1 ), aModelDim, aNumberOfElementsPerDirection );
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 )
                        {
                            //Possible side sets of children bottom left
                            // Children level = Parent level + 1
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 0 );
                            // +2, because col(0) = Level, col(1) = Element, col(2:5) = SideSet(1:4)
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tVar++;
                        }
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 1 ) > 0 )
                        {
                            //Possible side sets of children bottom right
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 1 );
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tVar++;
                        }
                        if ( tSideSet( i, 2 + 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 )
                        {
                            //Possible side sets of children top left
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 2 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tVar++;
                        }
                        if ( tSideSet( i, 2 + 1 ) > 0 || tSideSet( i, 2 + 2 ) > 0 )
                        {
                            //Possible side sets of children top right
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 3 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tVar++;
                        }
                    }
                    tSideSet.resize( tVar, 6 );
                }
                else if ( aModelDim == 3 )
                {
                    tHelp = 0;
                    for ( uint level = 0; level <= aLevel; level++ )
                    {
                        // 1 Sideset is always divided in 4 side sets in the children (1 Face of the cube)
                        tHelp += pow( 4, level );
                    }
                    // Level, Element, 6 Side sets = 8 cols
                    tSideSet.resize( Length * tHelp, 8 );
                    for ( uint i = 0; i < Length * ( tHelp - pow( 4, aLevel ) ); i++ )
                    {
                        if ( tSideSet( i, 0 ) >= aLevel )
                        {
                            //If finest level of updated side sets is reached, stop that loop
                            break;
                        }
                        tChildren = mBaseElement.give_children_of_element( tSideSet( i, 1 ), aModelDim, aNumberOfElementsPerDirection );
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 || tSideSet( i, 2 + 4 ) > 0 )
                        {
                            //Possible side sets of children bottom left in front (z-axis points in x-y plane)
                            // Children level = Parent level + 1
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 0 );
                            // +2, because col(0) = Level, col(1) = Element, col(2:7) = SideSet(1:6)
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tSideSet( tVar, 2 + 4 ) = tSideSet( i, 2 + 4 );
                            tSideSet( tVar, 2 + 5 ) = 0;
                            tVar++;
                        }
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 1 ) > 0 || tSideSet( i, 2 + 4 ) > 0 )
                        {
                            //Possible side sets of children bottom right in front (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 1 );
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tSideSet( tVar, 2 + 4 ) = tSideSet( i, 2 + 4 );
                            tSideSet( tVar, 2 + 5 ) = 0;
                            tVar++;
                        }
                        if ( tSideSet( i, 2 + 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 || tSideSet( i, 2 + 4 ) > 0)
                        {
                            //Possible side sets of children top left in front (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 2 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tSideSet( tVar, 2 + 4 ) = tSideSet( i, 2 + 4 );
                            tSideSet( tVar, 2 + 5 ) = 0;
                            tVar++;
                        }
                        if ( tSideSet( i, 2 + 1 ) > 0 || tSideSet( i, 2 + 2 ) > 0 || tSideSet( i, 2 + 4 ) > 0 )
                        {
                            //Possible side sets of children top right in front (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 3 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tSideSet( tVar, 2 + 4 ) = tSideSet( i, 2 + 4 );
                            tSideSet( tVar, 2 + 5 ) = 0;
                            tVar++;
                        }
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 || tSideSet( i, 2 + 5 ) > 0 )
                        {
                            //Possible side sets of children bottom left in back (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 4 );
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tSideSet( tVar, 2 + 4 ) = 0;
                            tSideSet( tVar, 2 + 5 ) = tSideSet( i, 2 + 5 );
                            tVar++;
                        }
                        if ( tSideSet( i, 2 ) > 0 || tSideSet( i, 2 + 1 ) > 0 || tSideSet( i, 2 + 5 ) > 0 )
                        {
                            //Possible side sets of children bottom right in back (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 5 );
                            tSideSet( tVar, 2 + 0 ) = tSideSet( i, 2 );
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = 0;
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tSideSet( tVar, 2 + 4 ) = 0;
                            tSideSet( tVar, 2 + 5 ) = tSideSet( i, 2 + 5 );
                            tVar++;
                        }

                        if ( tSideSet( i, 2 + 2 ) > 0 || tSideSet( i, 2 + 3 ) > 0 || tSideSet( i, 2 + 5 ) > 0 )
                        {
                            //Possible side sets of children top left in back (z-axis points in x-y plane)
                            tSideSet( tVar, 0 )     = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1 )     = tChildren( 6 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = 0;
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = tSideSet( i, 2 + 3 );
                            tSideSet( tVar, 2 + 4 ) = 0;
                            tSideSet( tVar, 2 + 5 ) = tSideSet( i, 2 + 5 );
                            tVar++;
                        }
                        if ( tSideSet( i, 2 + 1 ) > 0 || tSideSet( i, 2 + 2 ) > 0 || tSideSet( i, 2 + 5 ) > 0 )
                        {
                            //Possible side sets of children top right in back (z-axis points in x-y plane)
                            tSideSet( tVar, 0)      = tSideSet( i, 0 ) + 1;
                            tSideSet( tVar, 1)      = tChildren( 7 );
                            tSideSet( tVar, 2 + 0 ) = 0;
                            tSideSet( tVar, 2 + 1 ) = tSideSet( i, 2 + 1 );
                            tSideSet( tVar, 2 + 2 ) = tSideSet( i, 2 + 2 );
                            tSideSet( tVar, 2 + 3 ) = 0;
                            tSideSet( tVar, 2 + 4 ) = 0;
                            tSideSet( tVar, 2 + 5 ) = tSideSet( i, 2 + 5 );
                            tVar++;
                        }
                    }
                    tSideSet.resize( tVar, 8 );
                }
            }
        }
        Mat<uint> tMatSideSet;
        for ( uint i = 0; i < aNumSets( 1 ); i++ )
        {
            uint tVar = 0;
            if ( tSideSet.n_rows() > 0 )
            {
                tMatSideSet.set_size( tSideSet.n_rows(), 2 );
                for ( uint j = 0; j < tSideSet.n_rows(); j++ )
                {
                    // Check only side sets, which are in the range of the highest refinement level
                    if ( tSideSet(j,0) > aLevel )
                    {
                        break;
                    }
                    if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 2 ) == i + 1 )
                    {
                        tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                        tMatSideSet( tVar, 1 ) = 0;
                        tVar++;
                    }
                    if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 3 ) == i + 1 )
                    {
                        tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                        tMatSideSet( tVar, 1 ) = 1;
                        tVar++;
                    }
                    if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 4 ) == i + 1 )
                    {
                        tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                        tMatSideSet( tVar, 1 ) = 2;
                        tVar++;
                    }
                    if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 5 ) == i + 1 )
                    {
                        tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                        tMatSideSet( tVar, 1 ) = 3;
                        tVar++;
                    }
                    if ( aModelDim == 3 )
                    {
                        if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 6 ) == i + 1 )
                        {
                            tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                            tMatSideSet( tVar, 1 ) = 4;
                            tVar++;
                        }
                        if ( aElementActive.test( tSideSet( j, 1 ) ) == 1 && tSideSet( j, 7 ) == i + 1 )
                        {
                            tMatSideSet( tVar, 0 ) = tSideSet( j, 1 );
                            tMatSideSet( tVar, 1 ) = 5;
                            tVar++;
                        }
                    }
                }
            }
            if ( tVar > 0 )
            {
                tMatSideSet.resize( tVar, 2 );
                tSetsEntIds( tVarCell ) = tMatSideSet;
                tVarCell++;
            }
            else
            {
                tVarCell++;
            }
        }
    }
    return tSetsEntIds;
}

Cell<Mat<uint>>
Hierarchical_Mesh_Sets::update_blockset(
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aNumSets,
        Mat<uint> const & aBlockSet,
        Mat<uint> const & aElementListOnProc)
{
    // Initilize Cell for mesh output
    Cell<Mat<uint>> tSetsEntIds( aNumSets( 2 ) );
    Mat<uint> tBlockSetOutput( aElementListOnProc.length(), 1, 0 );
    uint tElementLevel = 0;
    uint tElementId = 0;
    //Update block set
    if ( aNumSets( 2 ) > 0 )
    {
        if ( aBlockSet.n_rows() > 0 )
        {
            for ( uint i = 0; i < aElementListOnProc.length(); i++ )
            {
                tElementId = aElementListOnProc( i );
                tElementLevel = mBaseElement.give_element_level( aElementListOnProc( i ), aModelDim, aNumberOfElementsPerDirection );
                if ( tElementLevel > 0 )
                {
                    //Check parent on level zero
                    tElementId = mBaseElement.give_parent_of_level_x( aElementListOnProc( i ), aModelDim, aNumberOfElementsPerDirection, 0);
                }
                //Grab block set id from element on level zero
                tBlockSetOutput( i ) = aBlockSet( tElementId );
            }
            //Temporary vector to save only block sets for each block set id
            Mat<uint> tBlockSetDummy;
            //Temporary variable for loop
            uint tVar = 0;
            for (uint j = 0; j <= aNumSets( 2 ); j++ )
            {
                tBlockSetDummy.set_size( aElementListOnProc.length(), 1, 0 );
                tVar = 0;
                for(uint i = 0; i < aElementListOnProc.length(); i++ )
                {
                    if( tBlockSetOutput( i ) == aNumSets( j ) )
                    {
                        tBlockSetDummy( tVar ) = aElementListOnProc( i );
                        tVar++;
                    }
                }
                tBlockSetDummy.resize( tVar, 1 );
                tSetsEntIds( j ) = tBlockSetOutput;
            }
        }
    }
    return tSetsEntIds;
}

