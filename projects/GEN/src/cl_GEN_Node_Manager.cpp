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

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::Node_Manager( mtk::Mesh* aMesh )
    {
        this->reset_background_nodes( aMesh );
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager::~Node_Manager()
    {
        this->delete_all_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::reset_background_nodes( mtk::Mesh* aMesh )
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

    uint Node_Manager::get_number_of_background_nodes()
    {
        return mBackgroundNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node_Manager::get_total_number_of_nodes()
    {
        return mBackgroundNodes.size() + mDerivedNodes.size();
    }

    //--------------------------------------------------------------------------------------------------------------

    const Node& Node_Manager::get_node( uint aNodeIndex )
    {
        if ( this->is_background_node( aNodeIndex ) )
        {
            return this->get_background_node( aNodeIndex );
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

    void Node_Manager::create_derived_node(
            const Vector< Background_Node* >& aBackgroundNodes,
            const Matrix< DDRMat >&           aParametricCoordinates,
            mtk::Geometry_Type                aGeometryType,
            mtk::Interpolation_Order          aInterpolationOrder )
    {
        mDerivedNodes.push_back(
                new Derived_Node(
                        this->get_total_number_of_nodes(),
                        aBackgroundNodes,
                        aParametricCoordinates,
                        aGeometryType,
                        aInterpolationOrder ) );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::add_derived_node( Derived_Node* aDerivedNode )
    {
        MORIS_ASSERT( aDerivedNode->get_index() == this->get_total_number_of_nodes(),
                "Attempted to add a derived node with node index %d when the next index should be %d",
                aDerivedNode->get_index(),
                this->get_total_number_of_nodes() );
        mDerivedNodes.push_back( aDerivedNode );
    }

    //--------------------------------------------------------------------------------------------------------------

    const Derived_Node& Node_Manager::get_derived_node( uint aDerivedNodeIndex ) const
    {
        MORIS_ASSERT( aDerivedNodeIndex >= mBackgroundNodes.size(),
                "A derived node was requested from the GEN node manager, but the index provided corresponds to a background node." );
        return *mDerivedNodes( aDerivedNodeIndex - mBackgroundNodes.size() );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::update_derived_node(
            uint        aNodeIndex,
            moris_id    aNodeID,
            moris_index aNodeOwner )
    {
        // Get derived node
        Derived_Node& tDerivedNode = this->get_derived_node( aNodeIndex );

        // Update node information
        tDerivedNode.set_id( aNodeID );
        tDerivedNode.set_owner( aNodeOwner );
    }

    //--------------------------------------------------------------------------------------------------------------

    real Node_Manager::get_node_coordinate_value(
            uint aNodeIndex,
            uint aCoordinateIndex )
    {
        return this->get_node( aNodeIndex ).get_coordinate_value( aCoordinateIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Node_Manager::node_depends_on_advs( uint aNodeIndex )
    {
        return this->get_node( aNodeIndex ).depends_on_advs();
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Node_Manager::get_number_of_derived_node_pdvs( uint aNodeIndex )
    {
        return this->get_derived_node( aNodeIndex ).get_num_pdvs();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::set_derived_node_starting_pdv_id( uint aNodeIndex, moris_id aStartingPDVID )
    {
        this->get_derived_node( aNodeIndex ).set_starting_pdv_id( aStartingPDVID );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id Node_Manager::get_derived_node_starting_pdv_id( uint aNodeIndex )
    {
        return this->get_derived_node( aNodeIndex ).get_starting_pdv_id();
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_id Node_Manager::get_derived_node_id( uint aNodeIndex )
    {
        return this->get_derived_node( aNodeIndex ).get_id();
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Node_Manager::get_derived_node_owner( uint aNodeIndex )
    {
        return this->get_derived_node( aNodeIndex ).get_owner();
    }

    //--------------------------------------------------------------------------------------------------------------

    void Node_Manager::append_dcoordinate_dadv_from_derived_node(
            uint                    aNodeIndex,
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor )
    {
        this->get_derived_node( aNodeIndex ).append_dcoordinate_dadv( aCoordinateSensitivities, aSensitivityFactor );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Node_Manager::get_coordinate_determining_adv_ids_from_derived_node( uint aNodeIndex )
    {
        return this->get_derived_node( aNodeIndex ).get_coordinate_determining_adv_ids();
    }

    //--------------------------------------------------------------------------------------------------------------

    Node_Manager& Node_Manager::get_trivial_instance()
    {
        static Node_Manager tManager( nullptr );
        return tManager;
    }

    //--------------------------------------------------------------------------------------------------------------

    Derived_Node& Node_Manager::get_derived_node( uint aDerivedNodeIndex )
    {
        MORIS_ASSERT( aDerivedNodeIndex >= mBackgroundNodes.size(),
                "A derived node was requested from the GEN node manager, but the index provided (%d) corresponds to a background node.",
                aDerivedNodeIndex );
        return *mDerivedNodes( aDerivedNodeIndex - mBackgroundNodes.size() );
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

}    // namespace moris::gen
