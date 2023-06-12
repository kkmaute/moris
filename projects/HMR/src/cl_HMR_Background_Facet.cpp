/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Background_Facet.cpp
 *
 */

#include "cl_HMR_Background_Facet.hpp"

#include "cl_HMR_Background_Element_Base.hpp"

namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------

        Background_Facet::Background_Facet( Background_Element_Base * aElementA,
                                            Background_Element_Base * aElementB,
                                            uint                      aIndexOnElementA )
        {
            if ( aElementA->is_padding() )
            {
                mLeaderElement = aElementB;
                mFollowerElement  = aElementA;
                mIndexOnLeader = this->get_index_on_other( aIndexOnElementA );
            }
            else if ( aElementB->is_padding() )
            {
                mLeaderElement = aElementA;
                mFollowerElement  = aElementB;
                mIndexOnLeader = aIndexOnElementA;
            }
            else if( aElementA->get_hmr_id() < aElementB->get_hmr_id() )
            {
                mLeaderElement = aElementA;
                mFollowerElement  = aElementB;
                mIndexOnLeader = aIndexOnElementA;
            }
            else
            {
                mLeaderElement = aElementB;
                mFollowerElement  = aElementA;
                mIndexOnLeader = this->get_index_on_other( aIndexOnElementA );
            }

            // set owner flags of leader and follower
        }

//-------------------------------------------------------------------------------

        Background_Facet::Background_Facet( Background_Element_Base * aElement,
                                            uint                      aIndexOnElement )
        {
            mLeaderElement = aElement;
            mFollowerElement  = nullptr;
            mIndexOnLeader = aIndexOnElement;
        }

//-------------------------------------------------------------------------------

        void Background_Facet::flag()
        {
            mFlag = true;
        }

//-------------------------------------------------------------------------------

        void Background_Facet::unflag()
        {
            mFlag = false;
        }

//-------------------------------------------------------------------------------

        bool Background_Facet::is_flagged() const
        {
            return mFlag;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base * Background_Facet::get_leader()
        {
            return mLeaderElement;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base * Background_Facet::get_follower()
        {
            return mFollowerElement;
        }

//-------------------------------------------------------------------------------

        void Background_Facet::set_follower( Background_Element_Base * aElement )
        {
            mFollowerElement = aElement;
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_leader() const
        {
            return mIndexOnLeader;
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_other( uint aIndex ) const
        {
            uint tNeighborFace[ 6 ] = { 2, 3, 0, 1, 5, 4 } ;
            return tNeighborFace[ aIndex ];
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_follower() const
        {
            return this->get_index_on_other( mIndexOnLeader );
        }

//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

