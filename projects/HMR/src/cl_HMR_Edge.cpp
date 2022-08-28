/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Edge.cpp
 *
 */

#include "cl_HMR_Edge.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace hmr
    {
        //------------------------------------------------------------------------------

        Edge::Edge( Mesh_Base       * aMesh,
                Background_Edge * aBackgroundEdge )
        {
            this->find_master( aMesh, aBackgroundEdge );
        }

        //------------------------------------------------------------------------------

        moris_id Edge::get_id() const
        {
            MORIS_ASSERT( mID != gNoID,"Edge ID not initialized");

            return mID;
        }

        //------------------------------------------------------------------------------

        moris_id Edge::get_index() const
        {
            MORIS_ASSERT( mID != gNoIndex,"Edge index not initialized");

            return mIndex;
        }

        //------------------------------------------------------------------------------

        void Edge::set_index( const moris_index & aIndex )
        {
            mIndex = aIndex;
        }
        //------------------------------------------------------------------------------

        void Edge::set_id( const moris_id & aID )
        {
            mID = aID;
        }

        //------------------------------------------------------------------------------

        moris_id Edge::get_owner() const
        {
            MORIS_ASSERT( mOwner != gNoID,"Edge ownership not initialized");

            return mOwner;
        }

        //------------------------------------------------------------------------------

        void Edge::set_owner( const moris_id & aOwner )
        {
            mOwner = aOwner;
        }

        //------------------------------------------------------------------------------

        moris::Cell< mtk::Vertex * > Edge::get_vertex_pointers() const
        {
            return moris::Cell< mtk::Vertex* >(0);
        }

        //------------------------------------------------------------------------------

        Matrix< IdMat > Edge::get_vertex_ids() const
        {
            // get number of vertices
            uint tNumberOfVertices =  this->get_number_of_vertices();

            // allocate output matrix
            Matrix< IdMat > aIDs( tNumberOfVertices, 1 );

            // loop over all vertices
            for( uint k=0; k<tNumberOfVertices; ++k )
            {
                // copy id into matrix
                aIDs( k ) = this->get_basis( k )->get_id();
            }

            // return ids
            return aIDs;
        }

        //------------------------------------------------------------------------------

        Matrix< IndexMat > Edge::get_vertex_inds() const
        {
            // get number of vertices
            uint tNumberOfVertices = this->get_number_of_vertices();

            // allocate output matrix
            Matrix< IndexMat > aIndices( tNumberOfVertices, 1 );

            // loop over all vertices
            for( uint k=0; k<tNumberOfVertices; ++k )
            {
                // copy id into matrix
                aIndices( k ) = this->get_basis( k )->get_index();
            }

            // return indices
            return aIndices;
        }

        //------------------------------------------------------------------------------

        Matrix< DDRMat > Edge::get_vertex_coords() const
        {
            // get number of vertices
            uint tNumberOfVertices = this->get_number_of_vertices();

            // allocate output matrix
            Matrix< DDRMat > aCoords( tNumberOfVertices, 3 );

            // loop over all basis
            for( uint k=0; k<tNumberOfVertices; ++k )
            {
                // fixme: do this in one line
                Matrix< DDRMat > tNodeCoords  = this->get_basis( k )->get_coords();

                // copy coords from vertex
                aCoords.set_row( k, trans( tNodeCoords ) );
            }

            return aCoords;
        }

        //------------------------------------------------------------------------------

        mtk::Geometry_Type Edge::get_geometry_type() const
        {
            return mtk::Geometry_Type::LINE;
        }

        //------------------------------------------------------------------------------

        void Edge::find_master( Mesh_Base       * aMesh,
                Background_Edge * aBackgroundEdge )
        {
            // master is element with lowest id and active

            // get activation pattern
            uint tPattern = aMesh->get_activation_pattern();

            // get number of elements form background edge
            luint tMinID = MORIS_LUINT_MAX;

            // ask background edge for number of elements
            uint tNumberOfElements = aBackgroundEdge->get_number_of_elements();

            // count elements that are neither deactive nor padding
            uint tCount = 0;

            for( uint k=0; k<tNumberOfElements; ++k )
            {
                // get element
                Background_Element_Base * tElement = aBackgroundEdge->get_element( k );
                if( ! tElement->is_deactive( tPattern ) && ! tElement->is_padding() )
                {
                    ++tCount;
                }
            }

            // allocate memory for element cells

            // allocate cell for element pointers
            mElements.resize( tCount, nullptr );
            mIndicesInElements.set_size( tCount, 1 );

            // populate arrays
            tCount = 0;
            for( uint k=0; k<tNumberOfElements; ++k )
            {
                // get element
                Background_Element_Base * tElement = aBackgroundEdge->get_element( k );
                if( ! tElement->is_deactive( tPattern ) && ! tElement->is_padding() )
                {
                    mElements( tCount )
                                   = aMesh->get_element_by_memory_index(
                                           tElement->get_memory_index() );
                    mIndicesInElements( tCount++ )
                    = aBackgroundEdge->get_index_on_element( k );
                }
            }

            // find master
            uint tKmin = MORIS_UINT_MAX;
            for( uint k=0; k<tCount; ++k )
            {
                // test if element is active
                if ( mElements( k )->is_active() &&  mElements( k )->
                        get_background_element()->get_hmr_id() < tMinID )
                {
                    tMinID = mElements( k )->get_background_element()->get_hmr_id();
                    tKmin = k;
                }
            }

            // test if no active master was found
            if( tKmin > tCount )
            {
                // search for refined element with lowest id
                for( uint k=0; k<tCount; ++k )
                {
                    if ( mElements( k )->is_refined() &&
                            mElements( k )->get_background_element()->get_hmr_id() < tMinID )
                    {
                        tMinID = mElements( k )->get_background_element()->get_hmr_id();
                        tKmin = k;
                    }
                }
            }

            MORIS_ASSERT( tKmin < tCount, "something went wrong while determining edge master" );

            // remember master index
            mIndexOfMaster = tKmin;

            // set owner of this edge
            mOwner = mElements( mIndexOfMaster )->get_owner();
        }
        //------------------------------------------------------------------------------

        bool Edge::is_active() const
        {
            uint tNumberOfElements = mElements.size();

            bool aActive = false;

            for( uint k=0; k<tNumberOfElements; ++k )
            {
                if( mElements( k )->is_active() )
                {
                    aActive = true;
                    break;
                }
            }
            return aActive;
        }

        //------------------------------------------------------------------------------

        uint Edge::get_number_of_elements() const
        {
            return mElements.size();
        }

        //------------------------------------------------------------------------------

        Element * Edge::get_element( const uint & aIndex )
        {
            return mElements( aIndex );
        }

        //------------------------------------------------------------------------------

        uint Edge::get_index_on_element(  const uint & aIndex ) const
        {
            return mIndicesInElements( aIndex );
        }

        //------------------------------------------------------------------------------

        Element * Edge::get_hmr_master()
        {
            return mElements( mIndexOfMaster );
        }

        //------------------------------------------------------------------------------

        uint Edge::get_index_on_master() const
        {
            return  mIndicesInElements( mIndexOfMaster );
        }

        //------------------------------------------------------------------------------
    } /* namespace hmr */
} /* namespace moris */

