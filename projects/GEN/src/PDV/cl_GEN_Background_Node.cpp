/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Background_Node.cpp
 *
 */

#include "cl_GEN_Background_Node.hpp"
#include "cl_GEN_Basis_Node.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Background_Node::Background_Node(
            uint                    aIndex,
            const Matrix< DDRMat >& aCoordinates )
            : Node( aIndex )
            , mCoordinates( aCoordinates )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Background_Node::get_global_coordinates() const
    {
        return mCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Background_Node::get_coordinate_value( uint aCoordinateIndex ) const
    {
        return mCoordinates( aCoordinateIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Cell< Basis_Node >        Background_Node::mDummyLocatorNodes = {};
    const Cell< Basis_Node >& Background_Node::get_locator_nodes() const
    {
        return mDummyLocatorNodes;
    }

    //--------------------------------------------------------------------------------------------------------------
}
