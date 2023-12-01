/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Parent_Derived_Node.cpp
 *
 */

#include "cl_GEN_Parent_Derived_Node.hpp"
#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Parent_Derived_Node::Parent_Derived_Node( Derived_Node* aDerivedNode )
            : mDerivedNode( aDerivedNode )
    {
    }

    //--------------------------------------------------------------------------------------------------------------
    
    uint Parent_Derived_Node::get_index() const
    {
        return mDerivedNode->get_index();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Matrix< DDRMat >& Parent_Derived_Node::get_global_coordinates() const
    {
        return mDerivedNode->get_global_coordinates();
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    const Matrix< DDRMat >& Parent_Derived_Node::get_parametric_coordinates() const
    {
        return mDerivedNode->get_parametric_coordinates();
    }
    
    //--------------------------------------------------------------------------------------------------------------

    Node* Parent_Derived_Node::get_node() const
    {
        return mDerivedNode;
    }

    //--------------------------------------------------------------------------------------------------------------
    
}
