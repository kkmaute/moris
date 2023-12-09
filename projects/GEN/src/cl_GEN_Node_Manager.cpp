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

    class Trivial_Node_Manager : public Node_Manager
    {
      public:
        /**
         * Trivial node manager constructor
         */
        Trivial_Node_Manager()
                : Node_Manager( nullptr )
        {
        };

        /**
         * When the geometry engine node manager is not yet created, derived nodes are not known yet.
         * Thus, the trivial version must reflect that all nodes are base nodes
         *
         * @param aNodeIndex Node index
         * @return True
         */
        bool is_base_node( uint aNodeIndex ) override
        {
            return true;
        }

        /**
         * Throws an error if a base node is requested from the trivial node manager, since this indicates that something is wrong.
         *
         * @param aIndex Node index
         * @return Error
         */
        Base_Node* get_base_node( uint aIndex ) override
        {
            MORIS_ERROR( false, "A base node was requested from a trivial node manager. An object is using the incorrect node manager." );
            return nullptr;
        }
    };

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

    uint Node_Manager::get_total_number_of_nodes()
    {
        return mBaseNodes.size() + mDerivedNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    Node* Node_Manager::get_node( uint aNodeIndex )
    {
        if ( this->is_base_node( aNodeIndex ) )
        {
            return this->get_base_node( aNodeIndex );
        }
        else
        {
            return this->get_derived_node( aNodeIndex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Node_Manager::is_base_node( uint aNodeIndex )
    {
        return aNodeIndex < mBaseNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    Base_Node* Node_Manager::get_base_node( uint aBaseNodeIndex )
    {
        return mBaseNodes( aBaseNodeIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::add_derived_node( Derived_Node* aDerivedNode )
    {
        MORIS_ASSERT( aDerivedNode->get_index() == this->get_total_number_of_nodes(),
                "Attempted to add a derived node with node index %d when the next index should be %d", aDerivedNode->get_index(), this->get_total_number_of_nodes() );
        mDerivedNodes.push_back( aDerivedNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node* Node_Manager::get_derived_node( uint aDerivedNodeIndex )
    {
        MORIS_ASSERT( aDerivedNodeIndex >= mBaseNodes.size(),
                "A derived node was requested from the GEN node manager, but the index provided corresponds to a base node." );
        return mDerivedNodes( aDerivedNodeIndex - mBaseNodes.size() );
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager& Node_Manager::get_trivial_instance()
    {
        static Trivial_Node_Manager tManager;
        return tManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::delete_all_nodes()
    {
        // Delete base nodes
        for ( Base_Node* iBaseNode : mBaseNodes )
        {
            delete iBaseNode;
        }
        mBaseNodes.clear();

        // Delete derived nodes
        for ( Derived_Node* iDerivedNode : mDerivedNodes )
        {
            delete iDerivedNode;
        }
        mDerivedNodes.clear();
    }

    //--------------------------------------------------------------------------------------------------------------

}
