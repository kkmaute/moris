/*
 * cl_Hierarchical_Mesh_Refinement.cpp
 *
 *  Created on: Dec 21, 2017
 *      Author: gleim
 */

#include "cl_Hierarchical_Mesh_Refinement.hpp" // STK/src/Hierarchical
using namespace moris;

void
Hierarchical_Mesh_Refinement::hierarchical_element_refinement(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> & aDeactivateElement,
        BoostBitset & aElementActive,
        BoostBitset & aElementRefined)
{
    //Check for uniqueness of elements and sort them
    aDeactivateElement = unique(aDeactivateElement);
    //Check for max element, which needs to be deactivated
    uint tMaxElement = aDeactivateElement.max();
    //Determine the children of this max element
    Mat<uint> tChildren = mBaseElement.give_children_of_element(tMaxElement,aModelDim,aNumberOfElementsPerDirection);
    tMaxElement = tChildren.max();
    //Compute the level, which will be the finest (highest)
    uint tLevel = mBaseElement.give_element_level(tMaxElement,aModelDim,aNumberOfElementsPerDirection);
    //Determine the number of elements, which are possible
    uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);
    //Resize the bitset if necessary
    if( aElementActive.size() < tNumberElements )
    {
        aElementActive.resize( tNumberElements );
        aElementRefined.resize( tNumberElements );
    }
    if ( aModelDim == 1 )
    {
        for ( uint i = 0; i < aDeactivateElement.length(); i++ )
        {
            //Children of the deactivated element
            tChildren = mBaseElement.give_children_of_element( aDeactivateElement( i ), aModelDim, aNumberOfElementsPerDirection );
            //Activate the children
            aElementActive.set( tChildren( 0 ) );
            aElementActive.set( tChildren( 1 ) );
            //Deactivate the element
            aElementActive.reset( aDeactivateElement( i ) );
            //Make the deactivated elemente to be passive
            aElementRefined.set( aDeactivateElement( i ) );
        }
    }
    else if ( aModelDim == 2 )
    {
        for ( uint i = 0; i < aDeactivateElement.length(); i++ )
        {
            //Children of the deactivated element
            tChildren = mBaseElement.give_children_of_element( aDeactivateElement( i ), aModelDim, aNumberOfElementsPerDirection );
            //Activate the children
            aElementActive.set( tChildren( 0 ) );
            aElementActive.set( tChildren( 1 ) );
            aElementActive.set( tChildren( 2 ) );
            aElementActive.set( tChildren( 3 ) );
            //Deactivate the element
            aElementActive.reset( aDeactivateElement( i ) );
            //Make the deactivated elemente to be passive
            aElementRefined.set( aDeactivateElement( i ) );
        }
    }
    else if ( aModelDim == 3 )
    {
        for ( uint i = 0; i < aDeactivateElement.length(); i++ )
        {
            //Children of the deactivated element
            tChildren = mBaseElement.give_children_of_element( aDeactivateElement( i ), aModelDim, aNumberOfElementsPerDirection );
            //Activate the children
            aElementActive.set( tChildren( 0 ) );
            aElementActive.set( tChildren( 1 ) );
            aElementActive.set( tChildren( 2 ) );
            aElementActive.set( tChildren( 3 ) );
            aElementActive.set( tChildren( 4 ) );
            aElementActive.set( tChildren( 5 ) );
            aElementActive.set( tChildren( 6 ) );
            aElementActive.set( tChildren( 7 ) );
            //Deactivate the element
            aElementActive.reset( aDeactivateElement( i ) );
            //Make the deactivated elemente to be passive
            aElementRefined.set( aDeactivateElement( i ) );
        }
    }
}

/* @TODO update_passive_outer_layer: switching loop and if statement should save some time */
/* @TODO check that this algorithm is correct */
void
Hierarchical_Mesh_Refinement::update_passive_outer_layer(
        uint const & aModelDim,
        uint const & aPolynomial,
        uint const & aLevel,
        Mat<uint> const & aNumberOfElementsPerDirection,
        Mat<uint> const & aElementRefinedListOnProc,
        BoostBitset & aElementRefined)
{
    //Save passive list in new variable
    Mat<uint> tElementList = aElementRefinedListOnProc;
    //Temporary variable for loop
    Mat<uint> tChildren;
    for( uint j = 0; j < aLevel; j++ )
    {
        Mat<uint> tElementChildrenList;
        if ( aModelDim == 1 )
        {
            tElementChildrenList.set_size( tElementList.length() * 2, 1, 0);
            for ( uint i = 0; i < tElementList.length(); i++ )
            {
                //Children of the deactivated element
                tChildren = mBaseElement.give_children_of_element( tElementList( i ), aModelDim, aNumberOfElementsPerDirection );
                //Save children of element in vector
                tElementChildrenList.rows( i * 2, ( i + 1 ) * 2 - 1 ) = tChildren.rows( 0, 1 );
                //Make the children to be passive
                for ( uint k = 0; k < 2; k++ )
                {
                    aElementRefined.set( tChildren( k ) );
                }
            }
        }
        else if ( aModelDim == 2 )
        {
            tElementChildrenList.set_size( tElementList.length() * 4, 1, 0);
            for ( uint i = 0; i < tElementList.length(); i++ )
            {
                //Children of the deactivated element
                tChildren = mBaseElement.give_children_of_element( tElementList( i ), aModelDim, aNumberOfElementsPerDirection );
                //Save children of element in vector
                tElementChildrenList.rows( i * 4, ( i + 1 ) * 4 - 1 ) = tChildren.rows( 0, 3 );
                //Make the children to be passive
                for ( uint k = 0; k < 4; k++ )
                {
                    aElementRefined.set( tChildren( k ) );
                }
            }
        }
        else if ( aModelDim == 3 )
        {
            tElementChildrenList.set_size( tElementList.length() * 8, 1, 0);
            for ( uint i = 0; i < tElementList.length(); i++ )
            {
                //Children of the deactivated element
                tChildren = mBaseElement.give_children_of_element( tElementList( i ), aModelDim, aNumberOfElementsPerDirection );
                //Save children of element in vector
                tElementChildrenList.rows( i * 8, ( i + 1 ) * 8 - 1 ) = tChildren.rows( 0, 7 );
                //Make the children to be passive
                for ( uint k = 0; k < 8; k++ )
                {
                    aElementRefined.set( tChildren( k ) );
                }
            }
        }
        //Childrens get in the loop the new element list
        tElementList = tElementChildrenList;
    }
}

//--------------------------------------------------------------------------------

void Hierarchical_Mesh_Refinement::activate_basisfunction(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        BoostBitset & aBasisActive,
        BoostBitset const & aElementActive,
        Mat<uint> const & aElementListOnProc,
        Mat<uint> & aActiveBasisList)
{
    Mat<uint> tBasis_of_element(pow(aPolynomial+1,aModelDim),1);

    // possible: we don't know yet if the element belnging to each basis function is active
    Mat<uint> tPossibleBasisFunction(aElementListOnProc.length()*pow(aPolynomial+1,aModelDim),1);
    uint tVarb = 0; // Temporary variable for loop
    for(uint i=0; i< aElementListOnProc.length(); i++)
    {
        tBasis_of_element = mHMRElement.give_basis_of_element(aElementListOnProc(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        tPossibleBasisFunction.rows(tVarb*pow(aPolynomial+1,aModelDim),(tVarb+1)*pow(aPolynomial+1,aModelDim)-1) = tBasis_of_element.rows(0,pow(aPolynomial+1,aModelDim)-1);
        tVarb++;
    }
    tPossibleBasisFunction=unique(tPossibleBasisFunction);
    aActiveBasisList.set_size(tPossibleBasisFunction.length(),1,0);
    uint tMaxBasis = tPossibleBasisFunction.max();
    uint tLevel = mBasis.give_basis_level(tMaxBasis,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
    uint tBasisTotal = mBasis.give_number_of_basis(tLevel,aModelDim,aPolynomial,aNumberOfElementsPerDirection);
    if( aBasisActive.size() < tBasisTotal )
    {
        aBasisActive.resize(tBasisTotal);
    }
    aBasisActive.reset();
    uint tParent = 0;
    uint tBasis_level;
    uint tVar = 0; // Temporary variable for loop
    Mat<uint> tWhichElementsActive(pow(aPolynomial+1,aModelDim),1);
    Mat<uint> tElement_of_basis(pow(aPolynomial+1,aModelDim),1);
    Mat<uint> tBasisOfElement;
    Mat<uint> tListOfDeactivatedBasis(pow(aPolynomial+1,aModelDim)*pow(aPolynomial+1,aModelDim),1);
    Mat<uint> tCountsOfNumbers; // Provides the counts of the number of values
    uint tHelp; // Temporary variable for the loop
    uint tSwitch = 0;
    uint tFurtherLevel;
    uint tVarTemp = 0;
    uint tParentOfParent;
    uint tSwitchOfParent;
    uint tVarActiveBasis = 0;

    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        tBasis_level = mBasis.give_basis_level(tPossibleBasisFunction(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        if( tBasis_level == 0 ) // Check the level of the element i
        {
            tElement_of_basis = mBasis.give_element_of_basis(tPossibleBasisFunction(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = aElementActive.test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                aBasisActive.set( tPossibleBasisFunction( i ) );
                aActiveBasisList( tVarActiveBasis ) = tPossibleBasisFunction( i );
                tVarActiveBasis++;
            }
        }
        else
        {
            break;
        }
    }
    //Activate basis on higher levels
    for(uint i=0; i<tPossibleBasisFunction.length(); i++) // Update basis functions from the list of possible basis functions
    {
        tSwitch = 0;
        tFurtherLevel = 0;
        tBasis_level = mBasis.give_basis_level(tPossibleBasisFunction(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection);
        if( tBasis_level > 0) // Check the level of the element i
        {
            tElement_of_basis = mBasis.give_element_of_basis(tPossibleBasisFunction(i),aModelDim,aPolynomial,aNumberOfElementsPerDirection);
            for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
            {
                tWhichElementsActive(j) = aElementActive.test(tElement_of_basis(j));
            }
            if ( sum( tWhichElementsActive ) > 0 )
            {
                tListOfDeactivatedBasis = tListOfDeactivatedBasis*0;
                while( tSwitch < 1)
                {
                    tSwitchOfParent = 0;
                    tVar = 0; // set to zero for the loop below
                    for(uint j=0; j<tElement_of_basis.length(); j++) // Update basis functions on each level, except level 0
                    {
                        tParent = mBaseElement.give_parent_of_element(tElement_of_basis(j),aModelDim,aNumberOfElementsPerDirection); // Parent of element
                        for(uint k = 0; k<tFurtherLevel; k++)
                        {
                            tParent = mBaseElement.give_parent_of_element(tParent,aModelDim,aNumberOfElementsPerDirection); // Parent of element
                        }
                        tParentOfParent = mBaseElement.give_parent_of_element(tParent,aModelDim,aNumberOfElementsPerDirection); // Parent of parent element
                        if( tParent == UINT_MAX )
                        {
                            tSwitch = 1;
                            tVar = tVarTemp;
                            break;
                        }
                        tBasisOfElement = mHMRElement.give_basis_of_element(tParent,aModelDim,aPolynomial,aNumberOfElementsPerDirection); // Get Basis functions of parent element
                        for(uint k=0; k<pow(aPolynomial+1,aModelDim); k++) // Update basis functions on each level, except level 0
                        {
                            if( aBasisActive.test(tBasisOfElement(k)) == 0 && aElementActive.test(tParent) == 0 )
                            {
                                tListOfDeactivatedBasis(tVar) = tBasisOfElement(k);
                                tVar++;
                            }
                            else
                            {
                                tSwitch = 1;
                            }
                        }
                        if( tParentOfParent != UINT_MAX && aElementActive.test(tParentOfParent) == 1)
                        {
                            tSwitchOfParent = 1;
                        }
                    }
                    tFurtherLevel++;
                    tVarTemp = tVar;
                    if( tSwitchOfParent == 1)
                    {
                        tSwitch = 0;
                    }
                }
                if( tVar > 0 && tSwitch == 1)
                {
                    /* @TODO find out why this crashes */
                    //Mat<uint> ReducedListOfDeactivatedBasis = tListOfDeactivatedBasis.rows(0,tVar-1);
                    Mat<uint> ReducedListOfDeactivatedBasis(tVar, 1);
                    ReducedListOfDeactivatedBasis.rows(0,tVar-1) = tListOfDeactivatedBasis.rows(0,tVar-1);

                    if (unique(ReducedListOfDeactivatedBasis).length() == 1)
                    {
                        ReducedListOfDeactivatedBasis.resize(ReducedListOfDeactivatedBasis.length()+1,1);
                    }
                    tCountsOfNumbers = histc(ReducedListOfDeactivatedBasis,unique(ReducedListOfDeactivatedBasis)); // Provides the counts of the number of values
                    Mat<uint> tHelpMat = (tCountsOfNumbers == pow(aPolynomial+1,aModelDim));
                    tHelp = sum(tHelpMat);
                    if ( tHelp > 0) // If all Elements with the support of the current basis function  are active or passive and at least one element is active, then the basis function is active
                    {
                        aBasisActive.set(tPossibleBasisFunction(i));
                        aActiveBasisList(tVarActiveBasis) = tPossibleBasisFunction(i);
                        tVarActiveBasis++;
                    }
                }
            }
        }
    }
    aActiveBasisList.resize(tVarActiveBasis,1);
}

//--------------------------------------------------------------------------------

/* @TODO check that this algorithm is correct */
void
Hierarchical_Mesh_Refinement::activate_basisfunction_new(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aLevel,
        BoostBitset & aBasisActive,
        BoostBitset & aBasisRefined,
        BoostBitset const & aElementActive,
        BoostBitset const & aElementRefined,
        Mat<uint> const & aElementListOnProc,
        Mat<uint> const & aElementRefinedListOnProc,
        Mat<uint> & aActiveBasisList)
{
    Mat<uint> tBasisOfElement;
    Mat<uint> tElementOfBasis;
    uint tNumberBasisInElement = pow( aPolynomial + 1, aModelDim );
    //Number of possible basis functions
    uint tNumberPossibleBasis = aElementListOnProc.length() * tNumberBasisInElement;
    //Save all active basis functions in a vector
    aActiveBasisList.set_size( tNumberPossibleBasis, 1, 0);
    //Temporary variable for the loop
    uint tSumCheck = 0;
    uint tCheckActive = 0;
    uint tVar = 0;
    //Check for finest level the number of basis functions
    uint tNumberOfBasis = mBasis.give_number_of_basis( aLevel, aModelDim, aPolynomial, aNumberOfElementsPerDirection );
    aBasisActive.resize( tNumberOfBasis );
    aBasisActive.reset();
    aBasisRefined.resize( tNumberOfBasis );
    aBasisRefined.reset();
    //Check for possible basis functions
    Mat<uint> tPossibleBasisFunction( aElementRefinedListOnProc.length() * tNumberBasisInElement, 1, 0 );
    for(uint i=0; i< aElementRefinedListOnProc.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element( aElementRefinedListOnProc( i ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        tPossibleBasisFunction.rows( tVar * tNumberBasisInElement, ( tVar + 1 ) * tNumberBasisInElement - 1 )
                                        = tBasisOfElement.rows( 0, tNumberBasisInElement - 1 );
        tVar++;
    }
    tPossibleBasisFunction = unique( tPossibleBasisFunction );
    tVar = 0;
    //Check which basis functions are passive (It it makes also some active basis functions passive, they will be activated in the second loop)
    for ( uint j = 0; j < tPossibleBasisFunction.length(); j++ )
    {
        tElementOfBasis = mBasis.give_element_of_basis( tPossibleBasisFunction( j ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        for( uint k = 0; k < tElementOfBasis.length(); k++ )
        {
            if ( aElementRefined.test( tElementOfBasis( k ) ) == 1 )
            {
                //Sum all elements which are passive
                tSumCheck++;
            }
        }
        if( tSumCheck == tNumberBasisInElement )
        {
            //If all elements are passive, the basis function is passive
            aBasisRefined.set( tPossibleBasisFunction( j ) );
        }
        //Reset to zero to use this variable for the next basis function
        tSumCheck = 0;
    }
    //Check for active basis functions
    tPossibleBasisFunction.set_size( aElementListOnProc.length() * tNumberBasisInElement, 1, 0 );
    tVar = 0;
    for ( uint i = 0; i < aElementListOnProc.length(); i++ )
    {
        tBasisOfElement = mHMRElement.give_basis_of_element( aElementListOnProc( i ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        tPossibleBasisFunction.rows( tVar * tNumberBasisInElement, ( tVar + 1 ) * tNumberBasisInElement - 1 )
                                        = tBasisOfElement.rows( 0, tNumberBasisInElement - 1 );
        tVar++;
    }
    tPossibleBasisFunction=unique(tPossibleBasisFunction);
    tVar = 0;
    //Check which basis functions are active
    for ( uint j = 0; j < tPossibleBasisFunction.length(); j++ )
    {
        tElementOfBasis = mBasis.give_element_of_basis( tPossibleBasisFunction( j ), aModelDim, aPolynomial, aNumberOfElementsPerDirection );
        for( uint k = 0; k < tElementOfBasis.length(); k++ )
        {
            if ( aElementRefined.test( tElementOfBasis( k ) ) == 1 )
            {
                //Sum all elements which are passive
                tSumCheck++;
            }
            else if( aElementActive.test( tElementOfBasis( k ) ) == 1 )
            {
                //Sum all elements which are active
                tSumCheck++;
                //At least one element needs to be active
                tCheckActive = 1;
            }
        }
        if(tSumCheck == tNumberBasisInElement && tCheckActive == 1 )
        {
            //The basis is active if all elements are passive and at least one is active
            aBasisActive.set( tPossibleBasisFunction( j ) );
            aActiveBasisList( tVar ) = tPossibleBasisFunction( j );
            tVar++;
        }
        tSumCheck = 0;
    }
    aActiveBasisList.resize( tVar, 1 );
}

//--------------------------------------------------------------------------------

Mat<uint>
Hierarchical_Mesh_Refinement::give_active_elements(
        uint const & aModelDim,
        uint const & aPolynomial,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aLevel,
        Mat<uint> const & aElementListOnProcInit,
        BoostBitset const & aElementActive) const
{
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level<aLevel+1; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aModelDim),level);
    }
    Mat<uint> tPossibleElementListOnProc(aElementListOnProcInit.length()*tMaxNumElemLevel,1,UINT_MAX);
    tMaxNumElemLevel -= pow(pow(2,aModelDim),aLevel);
    tPossibleElementListOnProc.rows(0,aElementListOnProcInit.length()-1) = aElementListOnProcInit.rows(0,aElementListOnProcInit.length()-1);
    uint tVar = aElementListOnProcInit.length();
    Mat<uint> tChildren(pow(2,aModelDim),1);
    //Mat<uint> tChildrenDummy(pow(2,aModelDim),1);
    for(uint i = 0; i<aElementListOnProcInit.length()*tMaxNumElemLevel; i++)
    {
        tChildren = mBaseElement.give_children_of_element(tPossibleElementListOnProc(i),aModelDim,aNumberOfElementsPerDirection);
        tPossibleElementListOnProc.rows(tVar+i*pow(2,aModelDim),tVar+(i+1)*pow(2,aModelDim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    tPossibleElementListOnProc = unique(tPossibleElementListOnProc);
    //    if( aElementActive.size() < (tPossibleElementListOnProc.max()+1)  )
    //        aElementActive.resize(tPossibleElementListOnProc.max()+1);
    Mat<uint> tElementListOnProc(tPossibleElementListOnProc.length(),1);
    tVar = 0;
    for(uint i = 0; i<tPossibleElementListOnProc.length(); i++)
    {
        if( aElementActive.test(tPossibleElementListOnProc(i)) == 1)
        {
            tElementListOnProc(tVar) = tPossibleElementListOnProc(i);
            tVar++;
        }
    }
    tElementListOnProc.resize(tVar,1);
    tElementListOnProc = sort(tElementListOnProc);
    return tElementListOnProc;
}

//--------------------------------------------------------------------------------

void Hierarchical_Mesh_Refinement::give_deactivated_elements(
        uint & aModelDim,
        uint & aPolynomial,
        Mat<uint> & aNumberOfElementsPerDirection,
        Mat<uint> & aDeactivateElements,
        BoostBitset & aElementActive,
        Mat<uint> & aElementListOnProcInit)
{
    //    uint tVar = 0; // Temporary variable for loop
    uint tLevel = 0;
    Mat<uint> tNeighbour; //Neighbors from an element
    if ( isempty(aDeactivateElements) == 0)
    {
        uint tMaxElement = (aDeactivateElements).max();
        Mat<uint> tChildren = mBaseElement.give_children_of_element(tMaxElement,aModelDim,aNumberOfElementsPerDirection);
        tMaxElement = tChildren.max();
        tLevel = mBaseElement.give_element_level(tMaxElement,aModelDim,aNumberOfElementsPerDirection);
        uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);
        if( aElementActive.size() < tNumberElements )
        {
            aElementActive.resize(tNumberElements);
        }
        for(uint i=0; i<(aDeactivateElements).length(); i++)
        {
            aElementActive.reset((aDeactivateElements)(i)); // Deactivate elements
        }
        //        tVar = (aDeactivateElements).length();
        //        (aDeactivateElements).resize(tVar+aElementActive.count(),1);
    }

    // Broadcast the deactivated elements to all procs for the first refinement level (It is needed to find dependencies at the border between procs)
    mMPI.broadcast_bitset_logical_and(aElementActive);

    //    if ( isempty(aDeactivateElements) == 0)
    //    {
    //        Mat<uint> ElementListOnProc = give_active_elements(aModelDim,aPolynomial,aNumberOfElementsPerDirection,tLevel,aElementListOnProcInit,aElementActive);
    //        //Check for linar dependencies
    //        if( aModelDim == 1 && isempty(ElementListOnProc) == 0 )
    //        {
    //            for(uint i=0; i<ElementListOnProc.length(); i++)
    //            {
    //                if( aElementActive.test(ElementListOnProc(i)) == 1 )
    //                {
    //                    uint tBuffer = 1; // Only one layer of neighbour elements are needed
    //                    tNeighbour = mBaseElement.give_neighbor_of_element(ElementListOnProc(i),aModelDim,tBuffer,aNumberOfElementsPerDirection);
    //                    if( aElementActive.test(tNeighbour(0)) == 0 && aElementActive.test(tNeighbour(1)) == 0 && aElementActive.test(tNeighbour(2)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(1);
    //                        aElementActive.reset(tNeighbour(1));
    //                        tVar++;
    //                    }
    //                }
    //            }
    //        }
    //        else if( aModelDim == 2 && isempty(ElementListOnProc) == 0 )
    //        {
    //            for(uint i=0; i<ElementListOnProc.length(); i++)
    //            {
    //                if( aElementActive.test(ElementListOnProc(i)) == 1 )
    //                {
    //                    uint tBuffer = 1; // Only one layer of neighbour elements are needed
    //                    tNeighbour = mBaseElement.give_neighbor_of_element(ElementListOnProc(i),aModelDim,tBuffer,aNumberOfElementsPerDirection);
    //                    if( aElementActive.test(tNeighbour(3)) == 0 && aElementActive.test(tNeighbour(5)) == 0 && aElementActive.test(tNeighbour(4)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(4);
    //                        aElementActive.reset(tNeighbour(4));
    //                        tVar++;
    //                    }
    //                    if( aElementActive.test(tNeighbour(1)) == 0 && aElementActive.test(tNeighbour(7)) == 0 && aElementActive.test(tNeighbour(4)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(4);
    //                        aElementActive.reset(tNeighbour(4));
    //                        tVar++;
    //                    }
    //                }
    //            }
    //        }
    //        else if( aModelDim == 3 && isempty(ElementListOnProc) == 0)
    //        {
    //            for(uint i=0; i<ElementListOnProc.length(); i++)
    //            {
    //                if( aElementActive.test(ElementListOnProc(i)) == 1 )
    //                {
    //                    uint tBuffer = 1; // Only one layer of neighbour elements are needed
    //                    tNeighbour = mBaseElement.give_neighbor_of_element(ElementListOnProc(i),aModelDim,tBuffer,aNumberOfElementsPerDirection);
    //                    if( aElementActive.test(tNeighbour(12)) == 0 && aElementActive.test(tNeighbour(14)) == 0 && aElementActive.test(tNeighbour(13)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(13);
    //                        aElementActive.reset(tNeighbour(13));
    //                        tVar++;
    //                    }
    //                    if( aElementActive.test(tNeighbour(10)) == 0 && aElementActive.test(tNeighbour(16)) == 0 && aElementActive.test(tNeighbour(13)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(13);
    //                        aElementActive.reset(tNeighbour(13));
    //                        tVar++;
    //                    }
    //                    if( aElementActive.test(tNeighbour(4)) == 0 && aElementActive.test(tNeighbour(22)) == 0 && aElementActive.test(tNeighbour(13)) == 1)
    //                    {
    //                        (aDeactivateElements)(tVar) = tNeighbour(13);
    //                        aElementActive.reset(tNeighbour(13));
    //                        tVar++;
    //                    }
    //                }
    //            }
    //        }
    //    }
    //    if ( isempty(aDeactivateElements) == 0)
    //    {
    //        aDeactivateElements.resize(tVar,1);
    //    }
}

/* @TODO make input parameters const */
/* @TODO consider using a BoostBitset for tElementStencilDensity  */
Mat<uint>
Hierarchical_Mesh_Refinement::find_elements_in_stencil(
        uint & aModelDim,
        uint & aPolynomial,
        Mat<uint> & aNumberOfElementsPerDirection,
        BoostBitset & aElementActive,
        Mat<uint> & aElementListOnProc,
        map<uint, real> & aElementField,
        real & aFeatureLowerBound,
        uint & aFeatureResolution,
        const uint      & aBufferElements,
        const bool      & aAdaptiveBufferLayer,
        const bool      & aStaircasebuffer )
{
    Mat<uint> tElementStencil;
    Mat<uint> tElementStencilDensity;
    Mat<uint> tElementNeighbours;
    Mat<uint> tRefineElement(aElementListOnProc.length(),1,0);
    Mat<uint> tRefineElementLvL(aElementListOnProc.length(),1,0);
    uint tVar = 0;
    uint tSum = 0;
    Mat<uint> tBasis;
    Mat<uint> tChildList;
    Mat<uint> tChildren;

    if( aFeatureResolution > 1 )
    {

        // increase boost bitset by one level
        uint tLevel = mBaseElement.give_element_level(
                aElementActive.size(),
                aModelDim,
                aNumberOfElementsPerDirection );
        tLevel++;
        uint tNumberOfElements = mBaseElement.give_number_of_elements(
                tLevel,
                aModelDim,
                aNumberOfElementsPerDirection);

        aElementActive.resize( tNumberOfElements );

        // This is temporary Bitset is only needed for searching of elements, which need to be deactivated
        BoostBitset tElementActiveDummy = aElementActive;

        //MORIS_ASSERT( aElementField.length() >= aElementListOnProc.max(),"Element field for find element in stencil has the wrong size");
        for(uint e = 0; e < aElementListOnProc.length(); e++)
        {
            if( aElementField.find(aElementListOnProc(e)) > aFeatureLowerBound )
            {
                tElementStencil = mHMRElement.give_neighbor_stencil_of_element(aFeatureResolution,aElementListOnProc(e),aModelDim,aPolynomial,aNumberOfElementsPerDirection,aElementActive);
                tElementStencilDensity.set_size(tElementStencil.n_rows(),tElementStencil.n_cols(),0);
                for(uint i = 0; i < tElementStencil.n_rows(); i++)
                {
                    for(uint j = 0; j < tElementStencil.n_cols(); j++)
                    {
                        if( tElementStencil(i,j) > 0 && tElementStencil(i,j) < UINT_MAX )
                        {
                            if( aElementField.find(tElementStencil(i,j)) < aFeatureLowerBound )
                            {
                                tElementStencilDensity(i,j) = 0;
                            }
                            else
                            {
                                tElementStencilDensity(i,j) = 1;
                            }
                        }
                    }
                }
                tElementNeighbours.set_size(1,tElementStencil.n_cols());
                Mat<uint> tCount(tElementStencil.n_rows(),1,0);
                for(uint i = 0; i < tElementStencil.n_rows(); i++)
                {
                    tSum = 1; // Count number of Elements with a density higher then DensityLowerBound
                    tElementNeighbours.row(0) = tElementStencil.row(i);

                    for(uint j = aFeatureResolution+2; j < tElementStencil.n_cols(); j++)
                    {
                        if( tElementStencilDensity(i,j) == 1 || tElementStencil(i,j) == UINT_MAX )
                        {
                            tSum++;
                        }
                        else
                        {
                            break;
                        }
                    }
                    for(uint j = aFeatureResolution; j > 0; j--)
                    {
                        if( tElementStencilDensity(i,j) == 1 || tElementStencil(i,j) == UINT_MAX )
                        {
                            tSum++;
                        }
                        else
                        {
                            break;
                        }
                    }
                    tCount(i) = tSum;
                }
                if( tCount.min() < aFeatureResolution)
                {
                    tRefineElement(tVar) = aElementListOnProc(e);

                    tElementActiveDummy.reset(tRefineElement(tVar));
                    tChildren = mBaseElement.give_children_of_element(
                            tRefineElement(tVar),
                            aModelDim,
                            aNumberOfElementsPerDirection);
                    for(uint l = 0; l < tChildren.length(); l++)
                        tElementActiveDummy.set(tChildren(l));

                    tVar++;
                }
            }
        }





        std::fprintf( stdout, "Number of Elements found in stencil %u\n", tVar);
        tRefineElement.resize(tVar,1);
        tRefineElement = unique(tRefineElement);

        // Parent elements need also be on the refinement list
        if( tVar > 0)
        {
            this->add_parents_to_list(aModelDim,aNumberOfElementsPerDirection,tRefineElement,tElementActiveDummy);

            uint tMaxRefinementElement = tRefineElement.max();
            uint tLevel = mBaseElement.give_element_level(tMaxRefinementElement,aModelDim,aNumberOfElementsPerDirection) + 1;
            uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);

            // Check if deactivated elements are on a higher level, then the bitset needs to be increased
            if( tElementActiveDummy.size() < tNumberElements )
            {
                tElementActiveDummy.resize(tNumberElements);
            }

            if ( aBufferElements > 0 )
            {
                this->create_buffer_layer(
                        aModelDim,
                        aNumberOfElementsPerDirection,
                        tLevel,
                        aBufferElements,
                        aAdaptiveBufferLayer,
                        aStaircasebuffer,
                        tRefineElement,
                        tElementActiveDummy);
            }
        }
    }

    return tRefineElement;
}

void Hierarchical_Mesh_Refinement::find_elements_for_refinement(
        uint & aModelDim,
        uint & aPolynomialDesign,
        Mat<uint> & aNumberOfElementsPerDirection,
        Mat<uint> & aElementListLastStepDesign,
        Mat<uint> & aDeactivateElements,
        Mat<real> & aElementField,
        BoostBitset & aElementActive,
        BoostBitset & aBasisActive,
        Mat<real> & aNodalField,
        bool & aTruncatedBsplines,
        real & aDensityLowerBound,
        real & aDensityIncrementBound,
        real & aDensityUpperBound,
        uint & aMaxDesignLevelOfRefinement,
        bool & aLvlSetMethod,
        uint & aMaxLevelSetLevelOfRefinement,
        uint & aBufferElements,
        bool & aAdaptiveBufferLayer,
        bool & aStaircasebuffer)
{
    uint tMaxElement = aElementListLastStepDesign.max();
    Mat<uint> tBasisOfElement = mHMRElement.give_basis_of_element(tMaxElement,aModelDim,aPolynomialDesign,aNumberOfElementsPerDirection);
    uint tMaxBasis = tBasisOfElement.max();
    Mat<uint> tChildren;
    uint tVar = 0;
    uint tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(aElementListLastStepDesign.length()-1),aModelDim,aNumberOfElementsPerDirection);
    tLevel++; // Add one level for the density of the children
    uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);
    BoostBitset tElementActiveDummy = aElementActive; // This is temporary Bitset is only needed for searching of elements, which need to be deactivated

    if( tElementActiveDummy.size() < tNumberElements )
        tElementActiveDummy.resize(tNumberElements);

    if( aElementField.length() < tNumberElements )
        aElementField.resize(tNumberElements,1);

    if( aNodalField.length() < tMaxBasis)
        aNodalField.resize(tMaxBasis,1);

    Mat<real> tKnotVector(2*(aPolynomialDesign+1),1,0);
    for(uint i = 0; i < tKnotVector.length(); i++)
        tKnotVector(i) = -(real)aPolynomialDesign + (real)i;
    Mat<real> tXi(aModelDim,1,0.5);
    Mat<uint> tElement_flag(aModelDim,1,0);
    Mat<real> Bspline = Bspline::build_spline_uniform_nd(tXi,tKnotVector,aPolynomialDesign,tElement_flag,aModelDim);
    Mat<real> tNodalLvLSet(pow(aPolynomialDesign+1,aModelDim),1); // Level set field of an element
    Mat<uint> tRefineElement(2*aElementListLastStepDesign.length(),1,0);
    // Check for element refinement. Check elemental density threshold and, if enabled, intersected elements of a level set field
    for(uint i = 0; i<aElementListLastStepDesign.length(); i++)
    {
        //        std::cout << " aElementListLastStepDesign(i) " << aElementListLastStepDesign(i) << " aElementField(aElementListLastStepDesign(i)) " << aElementField(aElementListLastStepDesign(i)) << std::endl;

        // get levelset for elements of this node
        tBasisOfElement = mHMRElement.give_basis_of_element(aElementListLastStepDesign(i),aModelDim,aPolynomialDesign,aNumberOfElementsPerDirection);
        for(uint j = 0; j < tBasisOfElement.length(); j++)
            tNodalLvLSet(j) = aNodalField(tBasisOfElement(j));

        // check if level set is less then zero and alement field within interval for this level
        if( tNodalLvLSet.min() < 0 && aElementField(aElementListLastStepDesign(i)) >= (aDensityLowerBound+aDensityIncrementBound*tLevel) && aElementField(aElementListLastStepDesign(i)) <= aDensityUpperBound    )
        {
            // copy field data to children
            tChildren = mBaseElement.give_children_of_element(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection);
            for(uint j = 0; j < tChildren.length(); j++)
                aElementField(tChildren(j)) = aElementField(aElementListLastStepDesign(i));

            // check if level of element is below max design level, and if so, mark for refinement
            tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection);
            if( tLevel < aMaxDesignLevelOfRefinement)
            {
                tRefineElement(tVar) = aElementListLastStepDesign(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else if( tLevel > 0 )
            {
                tLevel = aMaxDesignLevelOfRefinement-1;
                tRefineElement(tVar) = mBaseElement.give_parent_of_level_x(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
        // If this is an intersected element, this element needs to be refined
        if( aLvlSetMethod == true && tNodalLvLSet.min() < 0 && tNodalLvLSet.max() > 0)
        {

            tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection);

            // if element is below maximum refinement level, deactivate and activate children
            // otherwise, deactivate parent and activate all siblings
            if( tLevel < aMaxLevelSetLevelOfRefinement)
            {
                tRefineElement(tVar) = aElementListLastStepDesign(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else
            {
                tLevel = aMaxLevelSetLevelOfRefinement-1;
                tRefineElement(tVar) = mBaseElement.give_parent_of_level_x(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
    }
    tRefineElement.resize(tVar,1);
    tRefineElement = unique(tRefineElement);

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        this->add_parents_to_list(aModelDim,aNumberOfElementsPerDirection,tRefineElement,tElementActiveDummy);
    }

    if( tVar > 0 )
    {
        uint tMaxRefinementElement = tRefineElement.max();
        tLevel = mBaseElement.give_element_level(tMaxRefinementElement,aModelDim,aNumberOfElementsPerDirection) + 1;
        tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);
        if( tElementActiveDummy.size() < tNumberElements ) // Check if deactivated elements are on a higher level, then the bitset needs to be increased
            tElementActiveDummy.resize(tNumberElements);
    }
    // Add a buffer layer of elements
    if( tVar > 0 && aBufferElements > 0)
    {
        this->create_buffer_layer(aModelDim,aNumberOfElementsPerDirection,tLevel,aBufferElements,aAdaptiveBufferLayer,aStaircasebuffer,tRefineElement,tElementActiveDummy);
    }
    aDeactivateElements = tRefineElement;
}

void
Hierarchical_Mesh_Refinement::refinement_criteria_nodal_field(
        uint            & aModelDim,
        uint            & aPolynomialDesign,
        Mat<uint>       & aNumberOfElementsPerDirection,
        Mat<uint>       & aElementListLastStepDesign,
        Mat<uint>       & aDeactivateElements,
        BoostBitset     & aElementActive,
        BoostBitset     & aBasisActive,
        map<uint, real> & aNodalField,
        const bool      & aTruncatedBsplines,
        const real      & aNodalLowerBound,
        const real      & aNodalUpperBound,
        const uint      & aMaxDesignLevelOfRefinement,
        const uint      & aBufferElements,
        const bool      & aAdaptiveBufferLayer,
        const bool      & aStaircasebuffer)
{
    uint tMaxElement = aElementListLastStepDesign.max();
    Mat<uint> tBasisOfElement = mHMRElement.give_basis_of_element(tMaxElement,aModelDim,aPolynomialDesign,aNumberOfElementsPerDirection);

    Mat<uint> tChildren;
    uint tVar = 0;

    uint tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(aElementListLastStepDesign.length()-1),aModelDim,aNumberOfElementsPerDirection);
    tLevel++; // Add one level for the density of the children
    uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);

    if( aElementActive.size() < tNumberElements )
        aElementActive.resize(tNumberElements);

    // This is temporary Bitset is only needed for searching of elements, which need to be deactivated
    BoostBitset tElementActiveDummy = aElementActive;

    Mat<real> tKnotVector(2*(aPolynomialDesign+1),1,0);
    for(uint i = 0; i < tKnotVector.length(); i++)
    {
        tKnotVector(i) = -(real)aPolynomialDesign + (real)i;
    }
    
    Mat<real> tXi(aModelDim,1,0.5);
    Mat<uint> tElement_flag(aModelDim,1,0);
    Mat<real> Bspline = Bspline::build_spline_uniform_nd(tXi,tKnotVector,aPolynomialDesign,tElement_flag,aModelDim);
    Mat<real> tNodalLvLSet(pow(aPolynomialDesign+1,aModelDim),1); // Level set field of an element
    Mat<uint> tRefineElement(2*aElementListLastStepDesign.length(),1,0);

    // Check for element refinement. Check elemental density threshold and, if enabled, intersected elements of a level set field
    for(uint i = 0; i<aElementListLastStepDesign.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element(
                aElementListLastStepDesign(i),
                aModelDim,
                aPolynomialDesign,
                aNumberOfElementsPerDirection);

        for(uint j = 0; j < tBasisOfElement.length(); j++)
        {
            tNodalLvLSet(j) = aNodalField.find(tBasisOfElement(j));
        }

        // If elements crosses window, it must be refined
        if( tNodalLvLSet.min() < aNodalLowerBound && tNodalLvLSet.max() > aNodalUpperBound)
        {
            tLevel = mBaseElement.give_element_level(
                    aElementListLastStepDesign(i),
                    aModelDim,
                    aNumberOfElementsPerDirection);

            // if element is below maximum refinement level, deactivate and activate children
            // otherwise, deactivate parent and activate all siblings
            if( tLevel < aMaxDesignLevelOfRefinement)
            {
                tRefineElement(tVar) = aElementListLastStepDesign(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else if( tLevel > 0 )
            {
                tLevel = aMaxDesignLevelOfRefinement-1;
                tRefineElement(tVar) = mBaseElement.give_parent_of_level_x(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                {
                    tElementActiveDummy.set(tChildren(l));
                }
                tVar++;
            }
        }
    }
    tRefineElement.resize(tVar,1);
    tRefineElement = unique(tRefineElement);

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        this->add_parents_to_list(aModelDim,aNumberOfElementsPerDirection,tRefineElement,tElementActiveDummy);
    }

    if( tVar > 0 )
    {
        uint tMaxRefinementElement = tRefineElement.max();

        tLevel = mBaseElement.give_element_level(tMaxRefinementElement,aModelDim,aNumberOfElementsPerDirection) + 1;

        tNumberElements = mBaseElement.give_number_of_elements(
                tLevel,
                aModelDim,
                aNumberOfElementsPerDirection);

        // Check if deactivated elements are on a higher level, then the bitset needs to be increased
        if( tElementActiveDummy.size() < tNumberElements )
        {
            tElementActiveDummy.resize(tNumberElements);
        }
    }
    // Add a buffer layer of elements
    if( tVar > 0 && aBufferElements > 0)
    {
        this->create_buffer_layer(aModelDim,aNumberOfElementsPerDirection,tLevel,aBufferElements,aAdaptiveBufferLayer,aStaircasebuffer,tRefineElement,tElementActiveDummy);
    }
    aDeactivateElements = tRefineElement;
}

/* void Hierarchical_Mesh_Refinement::finish_refinement_list(
        Mat< uint >           & aDeactivateElements,
        const BoostBitset     & aElementActive,
        const uint            & aModelDim,
        const uint            & aPolynomialDesign,
        const Mat<uint>       & aNumberOfElementsPerDirection,
        const uint            & aBufferElements,
        const bool            & aAdaptiveBufferLayer,
        const bool            & aStaircasebuffer)
{

    if ( aDeactivateElements.length() > 0 )
    {
        // backup active list
        BoostBitset tElementActiveDummy = aElementActive;

        // calculate next level
        uint tMaxLevel = mBaseElement.give_element_level(
                aElementActive.size()-1,
                aModelDim,
                aNumberOfElementsPerDirection);

        // matrix containing children to be activated
        Mat< uint > tChildren;

        // loop over all elements in list
        for( uint k=0; k<aDeactivateElements.length(); ++k )
        {

            // deactivate this element
            tElementActiveDummy.reset( aDeactivateElements( k ) );

            uint tLevel = mBaseElement.give_element_level(
                    aDeactivateElements( k ),
                    aModelDim,
                    aNumberOfElementsPerDirection);

            if ( tLevel <  tMaxLevel )
            {
                // activate children
                tChildren = mBaseElement.give_children_of_element(
                        aDeactivateElements( k ),
                        aModelDim,
                        aNumberOfElementsPerDirection);

                // loop over all children
                for(uint i = 0; i < tChildren.length(); ++i)
                {
                    // flag children for activation
                    tElementActiveDummy.set( tChildren(i) );
                }
            }
        }

        // add parents to refinement list
        this->add_parents_to_list(
                aModelDim,
                aNumberOfElementsPerDirection,
                aDeactivateElements,
                tElementActiveDummy);

        // get maximum ID of element to refine
        //uint tMaxRefinementElement = aDeactivateElements.max();

        // create buffer for this element
        if( aBufferElements > 0)
        {
            this->create_buffer_layer(
                    aModelDim,
                    aNumberOfElementsPerDirection,
                    tMaxLevel,
                    aBufferElements,
                    aAdaptiveBufferLayer,
                    aStaircasebuffer,
                    aDeactivateElements,
                    tElementActiveDummy);
        }
    }
} */

//-------------------------------------------------------------------------------

/* @TODO make input parameters const */
void Hierarchical_Mesh_Refinement::refinement_criteria_element_field(
        uint            & aModelDim,
        uint            & aPolynomialDesign,
        Mat<uint>       & aNumberOfElementsPerDirection,
        Mat<uint>       & aElementListLastStepDesign,
        Mat<uint>       & aDeactivateElements,
        map<uint, real> & aElementField,
        BoostBitset     & aElementActive,
        BoostBitset     & aBasisActive,
        map<uint, real> & aNodalField,
        const bool      & aTruncatedBsplines,
        const real      & aDensityLowerBound,
        const real      & aDensityIncrementBound,
        const real      & aDensityUpperBound,
        const uint      & aMaxDesignLevelOfRefinement,
        const uint      & aBufferElements,
        const bool      & aAdaptiveBufferLayer,
        const bool      & aStaircasebuffer)
{
    // Fixme: more comments are needed
    uint tMaxElement = aElementListLastStepDesign.max();
    Mat<uint> tBasisOfElement = mHMRElement.give_basis_of_element(tMaxElement,aModelDim,aPolynomialDesign,aNumberOfElementsPerDirection);
    //uint tMaxBasis = tBasisOfElement.max();
    Mat<uint> tChildren;
    uint tVar = 0;
    uint tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(aElementListLastStepDesign.length()-1),aModelDim,aNumberOfElementsPerDirection);
    tLevel++; // Add one level for the density of the children
    uint tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);

    if( aElementActive.size() < tNumberElements )
    {
        aElementActive.resize(tNumberElements);
    }

    // This is temporary Bitset is only needed for searching of elements, which need to be deactivated
    BoostBitset tElementActiveDummy = aElementActive;

    /*if( aElementField.length() < tNumberElements )
        aElementField.resize(tNumberElements,1);

    if( aNodalField.length() < tMaxBasis)
        aNodalField.resize(tMaxBasis,1); */
    // Natural coordinate in the middle of the element
    Mat<real> tXi(aModelDim,1,0.5);

    // Compute the Bspline for the middle point
    Mat<real> Bspline = Bspline::build_spline_uniform_nd( tXi, aPolynomialDesign, aModelDim );

    // Level set field of an element
    Mat<real> tNodalLvLSet(pow(aPolynomialDesign+1,aModelDim),1);
    Mat<uint> tRefineElement(2*aElementListLastStepDesign.length(),1,0);

    // Check for element refinement. Check elemental density threshold and, if enabled, intersected elements of a level set field
    for(uint i = 0; i<aElementListLastStepDesign.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element(aElementListLastStepDesign(i),aModelDim,aPolynomialDesign,aNumberOfElementsPerDirection);
        for(uint j = 0; j < tBasisOfElement.length(); j++)
        {
            tNodalLvLSet(j) = aNodalField.find(tBasisOfElement(j));
        }

        uint tThisElement = aElementListLastStepDesign(i);

        real tEleDens = aElementField.find(tThisElement);
        // if nodal level set is below zero and element field intersects defined window
        if( tNodalLvLSet.min() < 0 && tEleDens >= (aDensityLowerBound+aDensityIncrementBound*tLevel) && tEleDens <= aDensityUpperBound )
        {
            tChildren = mBaseElement.give_children_of_element(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection);

            for(uint j = 0; j < tChildren.length(); j++)
            {
                aElementField[tChildren(j)] = tEleDens;
            }
            tLevel = mBaseElement.give_element_level(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection);

            // if element is below maximum refinement level, deactivate and activate children
            // otherwise, deactivate parent and activate all siblings
            if( tLevel < aMaxDesignLevelOfRefinement)
            {
                tRefineElement(tVar) = aElementListLastStepDesign(i);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
            else if( tLevel > 0 )
            {
                tLevel = aMaxDesignLevelOfRefinement-1;
                tRefineElement(tVar) = mBaseElement.give_parent_of_level_x(aElementListLastStepDesign(i),aModelDim,aNumberOfElementsPerDirection,tLevel);
                tElementActiveDummy.reset(tRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(tRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    tElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
    }

    tRefineElement.resize(tVar,1);
    tRefineElement = unique(tRefineElement);

    // Parent elements need also be on the refinement list
    if( tVar > 0)
    {
        this->add_parents_to_list(aModelDim,aNumberOfElementsPerDirection,tRefineElement,tElementActiveDummy);
    }

    //tRefineElement.print("tRefineElement");
    if( tVar > 0 )
    {
        uint tMaxRefinementElement = tRefineElement.max();
        tLevel = mBaseElement.give_element_level(tMaxRefinementElement,aModelDim,aNumberOfElementsPerDirection) + 1;
        tNumberElements = mBaseElement.give_number_of_elements(tLevel,aModelDim,aNumberOfElementsPerDirection);

        // Check if deactivated elements are on a higher level, then the bitset needs to be increased
        if( tElementActiveDummy.size() < tNumberElements )
        {
            tElementActiveDummy.resize(tNumberElements);
        }
    }

    // Add a buffer layer of elements
    if( tVar > 0 && aBufferElements > 0)
    {
        this->create_buffer_layer(
                aModelDim,
                aNumberOfElementsPerDirection,
                tLevel,
                aBufferElements,
                aAdaptiveBufferLayer,
                aStaircasebuffer,
                tRefineElement,
                tElementActiveDummy);
    }

    //tRefineElement.print("tRefineElement");
    aDeactivateElements = tRefineElement;
}

//-------------------------------------------------------------------------------

void Hierarchical_Mesh_Refinement::find_elements_of_object(
        const uint &      aModelDim,
        const uint &      aPolynomial,
        const uint &      aMaxLevelOfRefinement,
        const Mat<uint> & aNumberOfElementsPerDirection,
        BoostBitset &     aElementActive,
        const uint &      aBufferElements,
        const uint &      aBufferLevel,
        bool_t &          aAdaptiveBufferLayer,
        bool_t &          aStaircasebuffer,
        const Mat<uint> & aCandidateElementsToDeactivate,
        Mat<uint>&        aElementsToDeactivate)
{
    // get the number of elements
    uint tNumberOfElements
    = mBaseElement.give_number_of_elements(aMaxLevelOfRefinement,
            aModelDim,
            aNumberOfElementsPerDirection);

    // Fixme: this is a bit dirty and should done somewhere else
    // adapt BoostBitset length if necessary
    if( aElementActive.size() < tNumberOfElements)
    {
        aElementActive.resize(tNumberOfElements);
    }

    // create an internal copy of the BoostBitset
    BoostBitset tElementActiveInternal = aElementActive;

    // loop over all candidate elements
    // initialize output matrix
    aElementsToDeactivate.set_size(aCandidateElementsToDeactivate.length(), 1, UINT_MAX);
    uint tElementCounter = 0;
    Mat<uint> tChildren;
    for(uint i = 0; i < aCandidateElementsToDeactivate.length(); ++i)
    {
        // get this element number
        uint tThisElement = aCandidateElementsToDeactivate(i);

        // get level of current element
        uint tLevelOfThisElement = mBaseElement.give_element_level(
                tThisElement,
                aModelDim,
                aNumberOfElementsPerDirection);

        // check if element needs to be refined
        if( tLevelOfThisElement < aMaxLevelOfRefinement)
        {
            // add element to output list
            aElementsToDeactivate(tElementCounter) = tThisElement;

            // unset active flag of element
            tElementActiveInternal.reset(tThisElement);

            //Activate the children
            tChildren = mBaseElement.give_children_of_element(tThisElement,aModelDim,aNumberOfElementsPerDirection);
            for(uint l = 0; l < tChildren.length(); l++)
            {
                tElementActiveInternal.set(tChildren(l));
            }

            // increment counter
            ++tElementCounter;
        }
        else if( tLevelOfThisElement > 0 )
        {   // level could be equal or greater than max level
            // get parent on next higher level
            uint tParentOfThisElement = mBaseElement.give_parent_of_level_x(
                    tThisElement,
                    aModelDim,
                    aNumberOfElementsPerDirection,
                    aMaxLevelOfRefinement-1);

            // add parent on text higher level to list
            aElementsToDeactivate(tElementCounter) = tParentOfThisElement;

            // unset active flag of parent
            tElementActiveInternal.reset(tParentOfThisElement);

            //Activate the children
            tChildren = mBaseElement.give_children_of_element(tParentOfThisElement,aModelDim,aNumberOfElementsPerDirection);
            for(uint l = 0; l < tChildren.length(); l++)
            {
                tElementActiveInternal.set(tChildren(l));
            }

            // increment counter
            ++tElementCounter;
        }
    } // end loop over all elements

    // now we resize the output list
    aElementsToDeactivate.resize(tElementCounter, 1);

    if (tElementCounter > 0)
    {
        // we also need to add the parents of each element
        this->add_parents_to_list(
                aModelDim,
                aNumberOfElementsPerDirection,
                aElementsToDeactivate,
                tElementActiveInternal);

    }

    // now me make sure again, that the internal BoostBitset as correct length:
    if (tElementCounter > 0)
    {
        // get element with biggest ID
        auto tElementWithHighestNumber = aElementsToDeactivate.max();

        // get the level of this element plus one
        auto tNextLevel = mBaseElement.give_element_level(
                tElementWithHighestNumber, aModelDim, aNumberOfElementsPerDirection) + 1;

        // maximum possible size of BoostBitset
        auto tMaxNumberOfElements = mBaseElement.give_number_of_elements(tNextLevel,
                aModelDim,
                aNumberOfElementsPerDirection);
        // check for size
        if(tElementActiveInternal.size() < tMaxNumberOfElements)
        {
            tElementActiveInternal.resize(tMaxNumberOfElements);
        }
    }

    // finally, we add a buffer layer of elements
    if( tElementCounter > 0 && aBufferElements > 0)
    {
        this->create_buffer_layer(
                aModelDim,
                aNumberOfElementsPerDirection,
                aBufferLevel,
                aBufferElements,
                aAdaptiveBufferLayer,
                aStaircasebuffer,
                aElementsToDeactivate,
                tElementActiveInternal);
    }
}

//-------------------------------------------------------------------------------

/* @TODO make aElementListOnProcInit const */
Mat<uint>
Hierarchical_Mesh_Refinement::initial_refinement(
        uint const & aModelDim,
        Mat<uint> const & aNumberOfElementsPerDirection,
        uint const & aInitialRefinement,
        Mat<uint> & aElementListOnProcInit)
{
    uint tMaxNumElemLevel = 0;
    for(uint level = 0; level < aInitialRefinement; level ++)
    {
        tMaxNumElemLevel += pow(pow(2,aModelDim),level);
    }
    Mat<uint> tDeactiveElements((aElementListOnProcInit).length()*tMaxNumElemLevel,1,UINT_MAX);
    tMaxNumElemLevel -= pow(pow(2,aModelDim),aInitialRefinement-1);
    tDeactiveElements.rows(0,aElementListOnProcInit.length()-1) = aElementListOnProcInit.rows(0,aElementListOnProcInit.length()-1);
    uint tVar = aElementListOnProcInit.length();
    Mat<uint> tChildren(pow(2,aModelDim),1);
    Mat<uint> tChildrenDummy(pow(2,aModelDim),1);
    for(uint i = 0; i<aElementListOnProcInit.length()*tMaxNumElemLevel; i++)
    {
        tChildren = mBaseElement.give_children_of_element(tDeactiveElements(i),aModelDim,aNumberOfElementsPerDirection);
        tDeactiveElements.rows(tVar+i*pow(2,aModelDim),tVar+(i+1)*pow(2,aModelDim)-1) = tChildren.rows(0,tChildren.length()-1);
    }
    tDeactiveElements = unique(tDeactiveElements);
    return tDeactiveElements;
}

/* @TODO make input and output variables consistend with hpp file */
/* @TODO calculate tLevel only when needed */
/* @TODO rename tNeighbour to American English */
void Hierarchical_Mesh_Refinement::create_buffer_layer(
        const uint      & aModelDim,
        const Mat<uint> & aNumberOfElementsPerDirection,
        const uint      & aLevel,
        const uint      & aBufferElements,
        const bool      & aAdaptiveBufferLayer,
        const bool      & aStaircasebuffer,
        Mat<uint>       & aRefineElement,
        BoostBitset     & aElementActiveDummy)
{

    // get size of activation bitset
    uint tNumberOfPossibleElements = aElementActiveDummy.size();

    // list of elements to refine
    BoostBitset tElementsToRefine( tNumberOfPossibleElements );

    uint tLevel;
    Mat<uint> tChildren;
    Mat<uint> tNeighbour;

    // Starting point to save new elements, which need to be deactivated
    uint tLengthRefinement =  aRefineElement.length();
    uint tBuffer = 0;

    // loop over all elements ib aRefineElement
    for(uint j = 0; j < tLengthRefinement; j++)
    {
        // flag this element
        tElementsToRefine.set( aRefineElement( j ) );

        tLevel = mBaseElement.give_element_level(
                aRefineElement( j ),
                aModelDim,
                aNumberOfElementsPerDirection);

        // calculate neighbor element buffer
        if( aAdaptiveBufferLayer == true)
        {
            if( aStaircasebuffer == true)
            {
                tBuffer = aBufferElements;
            }
            else
            {
                tBuffer = aBufferElements*pow(2,tLevel);
            }
        }
        else
        {
            tBuffer = aBufferElements;
        }

        // get neighbors with respect to calculated buffer
        tNeighbour = mBaseElement.give_neighbor_of_element(
                aRefineElement( j ),
                aModelDim,
                tBuffer,
                aNumberOfElementsPerDirection);
        for(uint k = 0; k<tNeighbour.length(); k++)
        {
            // check if neighbor is within range and active
            if( tNeighbour(k) < aElementActiveDummy.size() &&  aElementActiveDummy.test(tNeighbour(k)) == 1 )
            {
                // flag neigbor in cue bitset
                tElementsToRefine.set( tNeighbour(k) );

                // unflag element in active bitset
                aElementActiveDummy.reset(tNeighbour(k));

                // get children of element
                tChildren = mBaseElement.give_children_of_element(
                        tNeighbour(k),
                        aModelDim,
                        aNumberOfElementsPerDirection);

                // flag children as active
                for(uint l = 0; l < tChildren.length(); l++)
                {
                    aElementActiveDummy.set(tChildren(l));
                }
            }
        }
    }

    // calculate how many elements are in cue
    uint tNumberOfElementsToRefine = tElementsToRefine.count();

    // resize output matrix
    aRefineElement.set_size( tNumberOfElementsToRefine, 1) ;

    // counter for output array
    uint tCount = 0;

    // loop over all possible elements
    for ( uint k=0; k<  tNumberOfPossibleElements; ++ k )
    {
        // test if element is flagged
        if ( tElementsToRefine.test( k ) )
        {
            // add element to output matrix
            aRefineElement( tCount ) = k;

            // increment counter
            ++tCount;
        }
    }

    // Adds the parent elements, which also need to be refined
    this->add_parents_to_list(
            aModelDim,
            aNumberOfElementsPerDirection,
            aRefineElement,
            aElementActiveDummy);
}
/* @TODO make input variables constant */
void Hierarchical_Mesh_Refinement::add_parents_to_list(
        const uint      & aModelDim,
        const Mat<uint> & aNumberOfElementsPerDirection,
        Mat<uint>       & aRefineElement,
        BoostBitset     & aElementActiveDummy)
{
    Mat<uint> tChildren;
    uint tMaxRefinementElement = aRefineElement.max();
    uint tLevel = mBaseElement.give_element_level(tMaxRefinementElement,aModelDim,aNumberOfElementsPerDirection);
    uint tVar = aRefineElement.length();
    uint tLengthRefinement = aRefineElement.length();
    aRefineElement.resize(aRefineElement.length()*(tLevel+1),1);

    // only elements below level 0 have parents
    if( tLevel > 0)
    {
        for(uint i = 0; i < tLengthRefinement; i++ )
        {
            tLevel = mBaseElement.give_element_level(aRefineElement(i),aModelDim,aNumberOfElementsPerDirection);
            for(uint j = 0; j < tLevel; j++)
            {
                aRefineElement(tVar) = mBaseElement.give_parent_of_level_x(aRefineElement(i),aModelDim,aNumberOfElementsPerDirection,j);
                aElementActiveDummy.reset(aRefineElement(tVar));
                tChildren = mBaseElement.give_children_of_element(aRefineElement(tVar),aModelDim,aNumberOfElementsPerDirection);
                for(uint l = 0; l < tChildren.length(); l++)
                    aElementActiveDummy.set(tChildren(l));
                tVar++;
            }
        }
    }
    aRefineElement.resize(tVar,1);
    aRefineElement = unique(aRefineElement);
}

//---------------------------------------------------------------------------------------------------------------------

/* @TODO make input variables constant */
void Hierarchical_Mesh_Refinement::update_basis_values_on_design_elements(
        uint            & aModelDim,
        uint            & aPolynomialDesign,
        Mat<uint>       & aNumberOfElementsPerDirection,
        Mat<uint>       & aElemLocaltoGlobalLastStepDesign,
        bool            & aTruncatedBsplines,
        BoostBitset     & aDesignBSplineActiveLastStep,
        map<uint, real> & aNodalField,
        Cell<Mat<real>> & aTMatrixParentChild)
{
    Mat<uint> tBasisOfElement;
    Mat<real> tTMatrix;
    Mat<uint> tIdField;
    map<uint, real> tNodalFieldNew = aNodalField;

    // Loop over all active elements of the last step
    for(uint i = 0; i < aElemLocaltoGlobalLastStepDesign.length(); i++)
    {
        // Get Tmatrix
        if( aTruncatedBsplines == false)
        {
            mTMatrix.give_Tmatrix_and_IdField(
                    aElemLocaltoGlobalLastStepDesign(i),
                    aModelDim,
                    aPolynomialDesign,
                    aNumberOfElementsPerDirection,
                    aDesignBSplineActiveLastStep,
                    tTMatrix,
                    tIdField,
                    aTMatrixParentChild);
        }
        else
        {
            mTMatrix.give_Truncated_Tmatrix_and_IdField(
                    aElemLocaltoGlobalLastStepDesign(i),
                    aModelDim,
                    aPolynomialDesign,
                    aNumberOfElementsPerDirection,
                    aDesignBSplineActiveLastStep,
                    tTMatrix,
                    tIdField,
                    aTMatrixParentChild);
        }

        // Get all basis of element at level of element (whether active or not)
        tBasisOfElement = mHMRElement.give_basis_of_element(
                aElemLocaltoGlobalLastStepDesign(i),
                aModelDim,
                aPolynomialDesign,
                aNumberOfElementsPerDirection);

        // Compute values at basis of current element
        for(uint j = 0; j < tTMatrix.n_cols(); j++)
        {
            tNodalFieldNew[tBasisOfElement(j)] = 0.0;
            for(uint k = 0; k < tTMatrix.n_rows(); k++)
            {
                tNodalFieldNew[tBasisOfElement(j)] += tTMatrix(k,j) * aNodalField.find(tIdField(k));
            }
        }
    }

    // update nodal field such that it contains values for all basis used by active elements in last mesh
    aNodalField = tNodalFieldNew;
}

//---------------------------------------------------------------------------------------------------------------------
//
// Computes elemental value by evaluating nodally defined field at element center
//

/* @TODO make input variables constant */
void Hierarchical_Mesh_Refinement::compute_element_field_on_design_elements(
        uint &             aModelDim,
        uint &             aPolynomialDesign,
        Mat<uint> &        aNumberOfElementsPerDirection,
        Mat<uint> &        aElemLocaltoGlobalLastStepDesign,
        bool &             aTruncatedBsplines,
        map<uint, real> &  aNodalField,
        map<uint, real> &  aElementField)
{
    aElementField.clear();

    Mat<uint> tBasisOfElement;
    Mat<real> tKnotVector(2*(aPolynomialDesign+1),1,0);

    for(uint i = 0; i < tKnotVector.length(); i++)
    {
        tKnotVector(i) = -(real)aPolynomialDesign + (real)i;
    }

    Mat<real> tXi(aModelDim,1,0.5);
    Mat<uint> tElement_flag(aModelDim,1,0);

    Mat<real> Bspline = Bspline::build_spline_uniform_nd(
            tXi,
            tKnotVector,
            aPolynomialDesign,
            tElement_flag,
            aModelDim);

    for(uint i = 0; i < aElemLocaltoGlobalLastStepDesign.length(); i++)
    {
        tBasisOfElement = mHMRElement.give_basis_of_element(
                aElemLocaltoGlobalLastStepDesign(i),
                aModelDim,aPolynomialDesign,
                aNumberOfElementsPerDirection);

        aElementField[aElemLocaltoGlobalLastStepDesign(i)] = 0.0;

        for(uint j = 0; j < tBasisOfElement.length(); j++)
        {
            aElementField[aElemLocaltoGlobalLastStepDesign(i)]  += Bspline(j) * aNodalField.find(tBasisOfElement(j));
        }
    }
}
