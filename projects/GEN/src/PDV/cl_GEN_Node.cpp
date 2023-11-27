/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Node.cpp
 *
 */

#include "cl_GEN_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Node::Node( uint aIndex )
            : mIndex( aIndex )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node::get_index()
    {
        return mIndex;
    }

    //--------------------------------------------------------------------------------------------------------------
}
