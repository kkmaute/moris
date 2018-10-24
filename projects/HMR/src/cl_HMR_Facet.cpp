#include "HMR_Globals.hpp"
#include "cl_HMR_Facet.hpp"
#include "cl_HMR_Mesh_Base.hpp"
#include "cl_HMR_Element.hpp"

namespace moris
{
    namespace hmr
    {
// ----------------------------------------------------------------------------

        Facet::Facet(
                Mesh_Base * aMesh,
                Background_Facet * aBackgroundFacet ) :
                mFacet( aBackgroundFacet )
        {
            // set pointer to master element
            mMaster = aMesh->get_element_by_memory_index(
                     aBackgroundFacet->get_master()->get_memory_index() );

            mIndexOnMaster = aBackgroundFacet->get_index_on_master();

            // set pointer to slave if background slave exists
            if( aBackgroundFacet->get_slave() != NULL )
            {
                mSlave =  aMesh->get_element_by_memory_index(
                        aBackgroundFacet->get_slave()->get_memory_index() );

                // fixme: this is not clean
                if( mMaster->get_basis( 0 ) == NULL )
                {
                    this->swap_master_and_slave();
                }
                MORIS_ASSERT( mMaster->get_basis( 0 ) != NULL, "Tried to create a facet without nodes");
            }
            else
            {
                mSlave = nullptr;
            }
        }

// ----------------------------------------------------------------------------

        Facet::~Facet(){}

// ----------------------------------------------------------------------------

        moris_id
        Facet::get_id() const
        {
            return mID;
        }

// ----------------------------------------------------------------------------

        moris_index
        Facet::get_index() const
        {
            return mIndex;
        }

// ----------------------------------------------------------------------------

        moris_id
        Facet::get_owner() const
        {
            return mMaster->get_owner();
        }

// ----------------------------------------------------------------------------

        moris::Cell< mtk::Vertex* >
        Facet::get_vertex_pointers() const
        {
            MORIS_ERROR( false, "get_vertex_pointers() not implemented for facet");
            return moris::Cell< mtk::Vertex* >( 0 );
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
            for( uint k=0; k<tNumberOfVertices; ++k )
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
            for( uint k=0; k<tNumberOfVertices; ++k )
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
            if( mSlave == NULL )
            {
                return mMaster->is_active();
            }
            else
            {
                return mMaster->is_active() || mSlave->is_active();
            }
        }

// ----------------------------------------------------------------------------

        void
        Facet::set_id( const moris_id & aID )
        {
            mID = aID;
        }

// ----------------------------------------------------------------------------

        void
        Facet::set_index( const moris_index & aIndex )
        {
            mIndex = aIndex;
        }

// ----------------------------------------------------------------------------

        uint
        Facet::get_index_on_master() const
        {
            return mIndexOnMaster;
        }

// ----------------------------------------------------------------------------

        uint
        Facet::get_index_on_slave() const
        {
            return mFacet->get_index_on_other( mIndexOnMaster );
        }

// ----------------------------------------------------------------------------
//  private:
// ----------------------------------------------------------------------------

        uint
        Facet::get_level() const
        {
            return mMaster->get_level();
        }

// ----------------------------------------------------------------------------

        Element *
        Facet::get_hmr_master()
        {
            return mMaster;
        }

// ----------------------------------------------------------------------------

        Element *
        Facet::get_hmr_slave()
        {
            return mSlave;
        }

// ----------------------------------------------------------------------------
        mtk::Cell *
        Facet::get_master()
        {
            return mMaster;
        }

// ----------------------------------------------------------------------------

        const mtk::Cell *
        Facet::get_master() const
        {
            return mMaster;
        }

// ----------------------------------------------------------------------------

        mtk::Cell *
        Facet::get_slave()
        {
            return mSlave;
        }

// ----------------------------------------------------------------------------

        const mtk::Cell *
        Facet::get_slave() const
        {
            return mSlave;
        }

// ----------------------------------------------------------------------------

        void
        Facet::swap_master_and_slave()
        {
            mIndexOnMaster = mFacet->get_index_on_other( mIndexOnMaster );
            Element * tSwap = mMaster;
            mMaster = mSlave;
            mSlave = tSwap;

        }
// ----------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */
