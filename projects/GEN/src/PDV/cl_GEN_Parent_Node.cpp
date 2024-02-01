/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Node.cpp
 *
 */

#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Parent_Node::Parent_Node(
            const Node&             aNode,
            const Matrix< DDRMat >& aParametricCoordinates )
            : mNode( aNode )
            , mParametricCoordinates( aParametricCoordinates )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parent_Node::get_index() const
    {
        return mNode.get_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Parent_Node::get_global_coordinates() const
    {
        return mNode.get_global_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Parent_Node::get_parametric_coordinates() const
    {
        return mParametricCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

}
