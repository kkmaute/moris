#include "cl_HMR_Background_Facet.hpp"
#include "cl_HMR_Background_Element_Base.hpp"

namespace moris
{
    namespace hmr
    {
//-------------------------------------------------------------------------------

        Background_Facet::Background_Facet(
                Background_Element_Base * aElementA,
                Background_Element_Base * aElementB,
                const  uint             & aIndexOnElementA  )
        {

            if ( aElementA->is_padding() )
            {
                mMasterElement = aElementB;
                mSlaveElement  = aElementA;
                mIndexOnMaster = this->get_index_on_other( aIndexOnElementA );
                mOwner         = aElementB->get_owner();
            }
            else if ( aElementB->is_padding() )
            {
                mMasterElement = aElementA;
                mSlaveElement  = aElementB;
                mIndexOnMaster = aIndexOnElementA;
                mOwner         = aElementA->get_owner();
            }
            else if( aElementA->get_domain_id() < aElementB->get_domain_id() )
            {
                mMasterElement = aElementA;
                mSlaveElement  = aElementB;
                mIndexOnMaster = aIndexOnElementA;
                mOwner         = aElementA->get_owner();
            }
            else
            {
                mMasterElement = aElementB;
                mSlaveElement  = aElementA;
                mIndexOnMaster = this->get_index_on_other( aIndexOnElementA );
                mOwner         = aElementB->get_owner();
            }
        }

//-------------------------------------------------------------------------------

        Background_Facet::Background_Facet(
                        Background_Element_Base * aElement,
                        const  uint             & aIndexOnElement,
                        const  moris_id         & aProcID )
        {
            mMasterElement = aElement;
            mSlaveElement  = nullptr;
            mIndexOnMaster = aIndexOnElement;
            mOwner         = aProcID;
        }

//-------------------------------------------------------------------------------

        void
        Background_Facet::flag()
        {
            mFlag = true;
        }

//-------------------------------------------------------------------------------

        void
        Background_Facet::unflag()
        {
            mFlag = false;
        }

//-------------------------------------------------------------------------------

        bool
        Background_Facet::is_flagged() const
        {
            return mFlag;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base *
        Background_Facet::get_master()
        {
            return mMasterElement;
        }

//-------------------------------------------------------------------------------

        Background_Element_Base *
        Background_Facet::get_slave()
        {
            return mSlaveElement;
        }

//-------------------------------------------------------------------------------

        void
        Background_Facet::set_slave( Background_Element_Base * aElement )
        {
            mSlaveElement = aElement;
        }

//-------------------------------------------------------------------------------

        uint
        Background_Facet::get_index_on_master() const
        {
            return mIndexOnMaster;
        }

//-------------------------------------------------------------------------------

        uint
        Background_Facet::get_index_on_other( const uint & aIndex ) const
        {
            uint tNeighborFace[ 6 ] = { 2, 3, 0, 1, 5, 4 } ;
            return tNeighborFace[ aIndex ];
        }

//-------------------------------------------------------------------------------

        uint
        Background_Facet::get_index_on_slave() const
        {
            return this->get_index_on_other( mIndexOnMaster );
        }

//-------------------------------------------------------------------------------

        moris_id
        Background_Facet::get_owner() const
        {
            return mOwner;
        }

//-------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
