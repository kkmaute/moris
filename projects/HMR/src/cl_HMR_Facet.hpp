/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Facet.hpp
 *
 */

#ifndef SRC_HMR_CL_HMR_FACE_BASE_HPP_
#define SRC_HMR_CL_HMR_FACE_BASE_HPP_

#include "cl_HMR_Background_Facet.hpp"
#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Facet.hpp"

namespace moris::hmr
{
    class Mesh_Base;

    class Element;
    class Basis;
    // ----------------------------------------------------------------------------

    class Facet : public mtk::Facet
    {
        // ----------------------------------------------------------------------------

      protected:

        // ----------------------------------------------------------------------------

        const Background_Facet* mFacet;

        // index of this face
        moris_id mID = gNoID;

        // id of this face
        moris_index mIndex = gNoIndex;

        // pointer to leader element
        Element* mLeader = nullptr;

        // pointer to follower element
        Element* mFollower = nullptr;

        uint mIndexOnLeader = MORIS_UINT_MAX;

        // pointer to parent
        // Facet * mParent = nullptr;

        // flag telling if facet has children
        // bool mHasChildren = false;

        // cell containing children
        // Vector< Facet * > mChildren;

        // child index of this facet
        // uint mChildIndex;

        // ----------------------------------------------------------------------------

      public:

        // ----------------------------------------------------------------------------

        Facet( Mesh_Base*         aMesh,
                Background_Facet* aBackgroundFacet );

        // ----------------------------------------------------------------------------

        virtual ~Facet();

        // ----------------------------------------------------------------------------

        /**
         * returns the domain wide id of the cell
         *
         * @return moris_id ID
         */
        moris_id get_id() const;

        // ----------------------------------------------------------------------------

        /**
         * returns the local index of the cell
         *
         * @return moris_index ID
         */
        moris_index get_index() const;

        // ----------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this cell
         */
        virtual uint get_number_of_vertices() const = 0;

        // ----------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this cell
         */
        moris_id get_owner() const;

        // ----------------------------------------------------------------------------

        /**
         * fills a Vector with pointers to connected vertices
         */
        Vector< mtk::Vertex* > get_vertex_pointers() const;

        // ----------------------------------------------------------------------------

        void remove_vertex_pointer( moris_index aIndex );

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with IDs of connected vertices
         */
        Matrix< IdMat > get_vertex_ids() const;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with Indices of connected vertices
         */
        Matrix< IndexMat > get_vertex_inds() const;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        virtual Matrix< DDRMat > get_vertex_coords() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns an enum that defines the geometry type of the element
         */
        virtual mtk::Geometry_Type get_geometry_type() const = 0;

        //------------------------------------------------------------------------------

        /**
         * returns the order of the element
         */
        virtual mtk::Interpolation_Order get_interpolation_order() const = 0;

        //------------------------------------------------------------------------------

        mtk::Cell* get_leader();

        //------------------------------------------------------------------------------

        const mtk::Cell* get_leader() const;

        //------------------------------------------------------------------------------

        mtk::Cell* get_follower();

        //------------------------------------------------------------------------------

        const mtk::Cell* get_follower() const;

        //-----------------------------------------------------------------------------

        uint get_index_on_leader() const;

        //-----------------------------------------------------------------------------

        uint get_index_on_follower() const;

        // ----------------------------------------------------------------------------

        virtual const mtk::Vertex* get_vertex( uint aIndex ) const = 0;

        // ----------------------------------------------------------------------------
        //      HMR public:
        // ----------------------------------------------------------------------------

        /**
         * inverts the order of the nodes
         */
        void flip();

        // ----------------------------------------------------------------------------

        bool is_active() const;

        // ----------------------------------------------------------------------------

        void set_id( const moris_id& aID );

        // ----------------------------------------------------------------------------

        /** increase id by given increment
         *
         * @param aIncrement increment by which id is increased
         */
        void increment_id( const moris_id aIncrement = 1 );

        // ----------------------------------------------------------------------------

        void set_index( const moris_index& aIndex );

        // ----------------------------------------------------------------------------

        Element* get_hmr_leader();

        // ----------------------------------------------------------------------------

        Element* get_hmr_follower();

        // ----------------------------------------------------------------------------

        uint get_level() const;

        // ----------------------------------------------------------------------------

        virtual Basis* get_basis( uint aIndex ) = 0;

        // ----------------------------------------------------------------------------

        virtual const Basis* get_basis( uint aIndex ) const = 0;

        // ----------------------------------------------------------------------------

      private:

        // ----------------------------------------------------------------------------

        void swap_leader_and_follower();

        // ----------------------------------------------------------------------------

    };  // class Facet

    // ----------------------------------------------------------------------------

} /* namespace moris */

#endif /* SRC_HMR_CL_HMR_FACE_BASE_HPP_ */
