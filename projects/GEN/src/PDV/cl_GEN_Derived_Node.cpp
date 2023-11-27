/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Derived_Node.cpp
 *
 */

#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Locator::Locator(
            Node* aNode,
            real aBasis )
            : mNode( aNode )
            , mBasis( aBasis )
    {
        MORIS_ASSERT( aNode, "A GEN Locator must be created with a valid node." );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Locator::get_index()
    {
        return mNode->get_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Locator::get_coordinates()
    {
        return mNode->get_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Locator::get_basis()
    {
        return mBasis;
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node::Derived_Node(
            uint            aIndex,
            Cell< Locator > aLocators )
            : Node( aIndex )
            , mLocators( aLocators )
    {
        MORIS_ASSERT( mLocators.size() > 0, "A derived GEN node must have at least one locator." );
        // TODO check for basis values here

        // Get contribution from first locator
        mCoordinates = Matrix< DDRMat >( 1, mLocators( 0 ).get_coordinates().length(), 0.0 );

        // Add contributions from other locators
        for ( auto iLocator : mLocators )
        {
            mCoordinates += iLocator.get_coordinates() * iLocator.get_basis();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Derived_Node::get_coordinates()
    {
        return mCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Locator >& Derived_Node::get_locators()
    {
        return mLocators;
    }

    //--------------------------------------------------------------------------------------------------------------

}
