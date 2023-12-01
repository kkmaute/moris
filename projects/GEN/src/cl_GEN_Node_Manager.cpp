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
#include "cl_GEN_Intersection_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::Node_Manager( mtk::Mesh* aMesh )
            : mBaseNodes( 0 )
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
            mBaseNodes.resize( aMesh->get_num_nodes() );
            for ( uint iNodeIndex = 0; iNodeIndex < mBaseNodes.size(); iNodeIndex++ )
            {
                mBaseNodes( iNodeIndex ) = new Base_Node( iNodeIndex, aMesh->get_node_coordinate( iNodeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node_Manager::get_number_of_base_nodes()
    {
        return mBaseNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node_Manager::get_total_number_of_nodes()
    {
        return mBaseNodes.size() + mDerivedNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    Base_Node* Node_Manager::get_base_node( uint aIndex )
    {
        return mBaseNodes( aIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::add_derived_node( Derived_Node* aDerivedNode )
    {
        mDerivedNodes.push_back( aDerivedNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node* Node_Manager::get_derived_node( uint aIndex )
    {
        MORIS_ASSERT( aIndex >= mBaseNodes.size(),
                "A derived node was requested from the GEN node manager, but the index provided corresponds to a base node." );
        return mDerivedNodes( aIndex - mBaseNodes.size() );
    }
    
    //--------------------------------------------------------------------------------------------------------------
    
    void Node_Manager::add_intersection_node( Intersection_Node* aIntersectionNode )
    {
        // All intersection nodes are derived nodes
        this->add_derived_node( aIntersectionNode );

        // Add to intersection node list
        mIntersectionNodes.push_back( aIntersectionNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Cell< Intersection_Node* >& Node_Manager::get_intersection_nodes()
    {
        return mIntersectionNodes;
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager& Node_Manager::get_trivial_instance()
    {
        static Node_Manager tManager( nullptr );
        return tManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::delete_all_nodes()
    {
        for ( Node* iNode : mBaseNodes )
        {
            delete iNode;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}
