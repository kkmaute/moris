/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Node.cpp
 *
 */

#include "cl_GEN_Base_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Base_Node::Base_Node(
            uint                    aIndex,
            const Matrix< DDRMat >& aCoordinates )
            : Node( aIndex )
            , mCoordinates( aCoordinates )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Base_Node::get_coordinates()
    {
        return mCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------
}
