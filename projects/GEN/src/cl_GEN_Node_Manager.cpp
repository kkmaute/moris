/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Node_Manager.cpp
 *
 */

#include "cl_GEN_Node_Manager.hpp"
#include "cl_MTK_Mesh_Core.hpp"
#include "cl_GEN_Base_Node.hpp"
#include "cl_GEN_Derived_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::Node_Manager( mtk::Mesh* aMesh )
            : mNodes( 0 )
            , mNumberOfBaseNodes( 0 )
    {
        this->reset_base_nodes( aMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::~Node_Manager()
    {
        this->delete_all_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::reset_base_nodes( mtk::Mesh* aMesh )
    {
        // Delete all old nodes
        this->delete_all_nodes();

        // Create new base GEN nodes if a mesh is given
        if ( aMesh )
        {
            mNumberOfBaseNodes = aMesh->get_num_nodes();
            mNodes.resize( mNumberOfBaseNodes );
            for ( uint iNodeIndex = 0; iNodeIndex < mNumberOfBaseNodes; iNodeIndex++ )
            {
                mNodes( iNodeIndex ) = new Base_Node( iNodeIndex, aMesh->get_node_coordinate( iNodeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Node* Node_Manager::get_node( uint aIndex ) const
    {
        return mNodes( aIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::add_derived_node( Derived_Node* aDerivedNode )
    {
        mNodes.push_back( aDerivedNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node* Node_Manager::get_derived_node( uint aIndex ) const
    {
        MORIS_ERROR( aIndex > mNumberOfBaseNodes, "A derived node was requested, but the index corresponds to a base node." );
        return static_cast< Derived_Node* >( mNodes( aIndex ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::delete_all_nodes()
    {
        for ( Node* iNode : mNodes )
        {
            delete iNode;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
