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

        Background_Facet::Background_Facet(       Background_Element_Base * aElementA,
                                                  Background_Element_Base * aElementB,
                                            const uint                    & aIndexOnElementA )
        {
            if ( aElementA->is_padding() )
            {
                mMasterElement = aElementB;
                mSlaveElement  = aElementA;
                mIndexOnMaster = this->get_index_on_other( aIndexOnElementA );
            }
            else if ( aElementB->is_padding() )
            {
                mMasterElement = aElementA;
                mSlaveElement  = aElementB;
                mIndexOnMaster = aIndexOnElementA;
            }
            else if( aElementA->get_hmr_id() < aElementB->get_hmr_id() )
            {
                mMasterElement = aElementA;
                mSlaveElement  = aElementB;
                mIndexOnMaster = aIndexOnElementA;
            }
            else
            {
                mMasterElement = aElementB;
                mSlaveElement  = aElementA;
                mIndexOnMaster = this->get_index_on_other( aIndexOnElementA );
            }

            // set owner flags of master and slave
        }

//-------------------------------------------------------------------------------

        Background_Facet::Background_Facet(       Background_Element_Base * aElement,
                                            const uint                    & aIndexOnElement )
        {
            mMasterElement = aElement;
            mSlaveElement  = nullptr;
            mIndexOnMaster = aIndexOnElement;
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

        Background_Element_Base * Background_Facet::get_master()
        {
            return mMasterElement;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base * Background_Facet::get_slave()
        {
            return mSlaveElement;
        }

//-------------------------------------------------------------------------------

        void Background_Facet::set_slave( Background_Element_Base * aElement )
        {
            mSlaveElement = aElement;
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_master() const
        {
            return mIndexOnMaster;
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_other( const uint & aIndex ) const
        {
            uint tNeighborFace[ 6 ] = { 2, 3, 0, 1, 5, 4 } ;
            return tNeighborFace[ aIndex ];
        }

//-------------------------------------------------------------------------------

        uint Background_Facet::get_index_on_slave() const
        {
            return this->get_index_on_other( mIndexOnMaster );
        }

//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

