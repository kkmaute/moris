/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Facet.cpp
 *
 */

#include "cl_HMR_Facet.hpp"

#include "cl_HMR_Element.hpp"
#include "cl_HMR_Mesh_Base.hpp"
#include "HMR_Globals.hpp"

namespace moris::hmr
{
    // ----------------------------------------------------------------------------

    Facet::Facet(
            Mesh_Base*        aMesh,
            Background_Facet* aBackgroundFacet )
            : mFacet( aBackgroundFacet )
    {
        // set pointer to leader element
        MORIS_ASSERT( aBackgroundFacet->get_leader() != nullptr, "Background leader facet is nullptr" );

        mLeader = aMesh->get_element_by_memory_index( aBackgroundFacet->get_leader()->get_memory_index() );

        mIndexOnLeader = aBackgroundFacet->get_index_on_leader();

        // set pointer to follower if background follower exists
        if ( aBackgroundFacet->get_follower() != nullptr )
        {
            mFollower = aMesh->get_element_by_memory_index( aBackgroundFacet->get_follower()->get_memory_index() );

            // fixme: this is not clean
            if ( mLeader->get_basis( 0 ) == nullptr )
            {
                this->swap_leader_and_follower();
            }
            MORIS_ASSERT( mLeader->get_basis( 0 ) != nullptr, "Tried to create a facet without nodes" );
        }
        else
        {
            mFollower = nullptr;
        }
    }

    // ----------------------------------------------------------------------------

    Facet::~Facet(){}

    // ----------------------------------------------------------------------------

    moris_id
    Facet::get_id() const
    {
        MORIS_ASSERT( mID != gNoID, "Facet ID not initialized" );

        return mID;
    }

    // ----------------------------------------------------------------------------

    moris_index
    Facet::get_index() const
    {
        MORIS_ASSERT( mID != gNoIndex, "Facet index not initialized" );

        return mIndex;
    }

    // ----------------------------------------------------------------------------

    moris_id
    Facet::get_owner() const
    {
        MORIS_ASSERT( mLeader->get_owner() != gNoID, "Faceted ownership not initialized" );

        return mLeader->get_owner();
    }

    // ----------------------------------------------------------------------------

    moris::Vector< mtk::Vertex* >
    Facet::get_vertex_pointers() const
    {
        MORIS_ERROR( false, "get_vertex_pointers() not implemented for facet" );
        return moris::Vector< mtk::Vertex* >( 0 );
    }

    // ----------------------------------------------------------------------------

    // TODO MESHCLEANUP
    void
    Facet::remove_vertex_pointer( moris_index aIndex )
    {
        std::cout << "In HMR Facet" << std::endl;
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Facet::get_vertex_ids() const
    {
        // get number of vertices
        uint tNumberOfVertices = this->get_number_of_vertices();

        // allocate output vector
        Matrix< IdMat > aIDs( tNumberOfVertices, 1 );

        // loop over all vertices
        for ( uint k = 0; k < tNumberOfVertices; ++k )
        {
            aIDs( k ) = this->get_vertex( k )->get_id();
        }

        return aIDs;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Facet::get_vertex_inds() const
    {
        // get number of vertices
        uint tNumberOfVertices = this->get_number_of_vertices();

        // allocate output vector
        Matrix< IndexMat > aInds( tNumberOfVertices, 1 );

        // loop over all vertices
        for ( uint k = 0; k < tNumberOfVertices; ++k )
        {
            aInds( k ) = this->get_vertex( k )->get_index();
        }

        return aInds;
    }
    // ----------------------------------------------------------------------------
    //  HMR public:
    // ----------------------------------------------------------------------------

    bool
    Facet::is_active() const
    {
        if ( mFollower == nullptr )
        {
            return mLeader->is_active();
        }
        else
        {
            return mLeader->is_active() || mFollower->is_active();
        }
    }

    // ----------------------------------------------------------------------------

    void
    Facet::set_id( const moris_id& aID )
    {
        mID = aID;
    }

    // ----------------------------------------------------------------------------

    void
    Facet::increment_id( const moris_id aIncrement )
    {
        mID += aIncrement;

        MORIS_ASSERT( mID >= 0, "After increment has been applied Facet ID still negative." );
    }

    // ----------------------------------------------------------------------------

    void
    Facet::set_index( const moris_index& aIndex )
    {
        mIndex = aIndex;
    }

    // ----------------------------------------------------------------------------

    uint
    Facet::get_index_on_leader() const
    {
        return mIndexOnLeader;
    }

    // ----------------------------------------------------------------------------

    uint
    Facet::get_index_on_follower() const
    {
        return mFacet->get_index_on_other( mIndexOnLeader );
    }

    // ----------------------------------------------------------------------------
    //  private:
    // ----------------------------------------------------------------------------

    uint
    Facet::get_level() const
    {
        return mLeader->get_level();
    }

    // ----------------------------------------------------------------------------

    Element*
    Facet::get_hmr_leader()
    {
        return mLeader;
    }

    // ----------------------------------------------------------------------------

    Element*
    Facet::get_hmr_follower()
    {
        return mFollower;
    }

    // ----------------------------------------------------------------------------
    mtk::Cell*
    Facet::get_leader()
    {
        return mLeader;
    }

    // ----------------------------------------------------------------------------

    const mtk::Cell*
    Facet::get_leader() const
    {
        return mLeader;
    }

    // ----------------------------------------------------------------------------

    mtk::Cell*
    Facet::get_follower()
    {
        return mFollower;
    }

    // ----------------------------------------------------------------------------

    const mtk::Cell*
    Facet::get_follower() const
    {
        return mFollower;
    }

    // ----------------------------------------------------------------------------

    void
    Facet::swap_leader_and_follower()
    {
        mIndexOnLeader = mFacet->get_index_on_other( mIndexOnLeader );
        Element* tSwap = mLeader;
        mLeader        = mFollower;
        mFollower         = tSwap;
    }

    // ----------------------------------------------------------------------------

} /* namespace moris */
