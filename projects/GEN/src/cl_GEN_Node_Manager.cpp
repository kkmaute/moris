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
#include "cl_GEN_Basis_Node.hpp"

namespace moris::ge
{

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::Node_Manager( mtk::Mesh* aMesh )
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
            // Set that a mesh was given
            mMeshGiven = true;

            // Get number of nodes from the mesh
            uint tNumberOfBackgroundNodes = aMesh->get_num_nodes();

            // Reserve space
            mBackgroundNodes.reserve( tNumberOfBackgroundNodes );

            // Populate vector
            for ( uint iNodeIndex = 0; iNodeIndex < tNumberOfBackgroundNodes; iNodeIndex++ )
            {
                mBackgroundNodes.emplace_back( iNodeIndex, aMesh->get_node_coordinate( iNodeIndex ) );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node_Manager::get_total_number_of_nodes()
    {
        return mBackgroundNodes.size() + mDerivedNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    Node* Node_Manager::get_node( uint aNodeIndex )
    {
        if ( this->is_background_node( aNodeIndex ) )
        {
            return &( this->get_background_node( aNodeIndex ) );
        }
        else
        {
            return this->get_derived_node( aNodeIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Node_Manager::is_background_node( uint aNodeIndex ) const
    {
        return aNodeIndex < mBackgroundNodes.size() or not mMeshGiven;
    }

    //--------------------------------------------------------------------------------------------------------------

    Background_Node& Node_Manager::get_background_node( uint aBackgroundNodeIndex )
    {
        return mBackgroundNodes( aBackgroundNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::add_derived_node( Derived_Node* aDerivedNode )
    {
        MORIS_ASSERT( aDerivedNode->get_index() == this->get_total_number_of_nodes(),
                "Attempted to add a derived node with node index %d when the next index should be %d", aDerivedNode->get_index(), this->get_total_number_of_nodes() );
        mDerivedNodes.push_back( aDerivedNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node* Node_Manager::get_derived_node( uint aDerivedNodeIndex ) const
    {
        MORIS_ASSERT( aDerivedNodeIndex >= mBackgroundNodes.size(),
                "A derived node was requested from the GEN node manager, but the index provided corresponds to a background node." );
        return mDerivedNodes( aDerivedNodeIndex - mBackgroundNodes.size() );
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
        // Clear background nodes
        mBackgroundNodes.clear();

        // Delete derived nodes
        for ( Derived_Node* iDerivedNode : mDerivedNodes )
        {
            delete iDerivedNode;
        }
        mDerivedNodes.clear();
    }

    //--------------------------------------------------------------------------------------------------------------

}
