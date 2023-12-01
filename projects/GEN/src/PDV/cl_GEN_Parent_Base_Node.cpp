/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Base_Node.cpp
 *
 */

#include "cl_GEN_Parent_Base_Node.hpp"
#include "cl_GEN_Base_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Parent_Base_Node::Parent_Base_Node(
            Base_Node*              aBaseNode,
            const Matrix< DDRMat >& aParametricCoordinates )
            : mBaseNode( aBaseNode )
            , mParametricCoordinates( aParametricCoordinates )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Parent_Base_Node::get_index() const
    {
        return mBaseNode->get_index();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Parent_Base_Node::get_global_coordinates() const
    {
        return mBaseNode->get_global_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Parent_Base_Node::get_parametric_coordinates() const
    {
        return mParametricCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    Node* Parent_Base_Node::get_node() const
    {
        return mBaseNode;
    }

    //--------------------------------------------------------------------------------------------------------------

}
