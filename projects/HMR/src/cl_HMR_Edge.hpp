/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Edge.hpp
 *
 */

#pragma once

#include "cl_HMR_Element.hpp"
#include "cl_HMR_Mesh_Base.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Edge.hpp"
#include "cl_MTK_Vertex.hpp"

namespace moris::hmr
{

    class Mesh_Base;
    // ----------------------------------------------------------------------------

    class Edge : public mtk::Edge
    {
        // id of this edge
        moris_id mID = gNoID;

        // index of this edge
        moris_index mIndex = gNoIndex;

        // pointer with indices in elements
        Matrix< DDUMat > mIndicesInElements;

        // owner of this edge
        moris_id mOwner = gNoID;

        // ----------------------------------------------------------------------------

      protected:
        // ----------------------------------------------------------------------------

        // index on leader
        uint mIndexOfLeader;

        // pointer with elements
        Vector< Element* > mElements;

        // ----------------------------------------------------------------------------

      public:
        // ----------------------------------------------------------------------------

        /**
         * constructor
         */
        Edge(
                hmr::Mesh_Base*  aMesh,
                Background_Edge* aBackgroundEdge );

        // ----------------------------------------------------------------------------

        /**
         * trivial destructor
         */
        ~Edge() override{};

        // ----------------------------------------------------------------------------

        /**
         * returns the domain wide id of the edge
         *
         * @return moris_id ID
         */
        moris_id get_id() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns the local index of the edge
         *
         * @return moris_index ID
         */
        moris_index get_index() const override;

        // ----------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this edge
         */
        uint get_number_of_vertices() const override = 0;

        // ----------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this edge
         * ( this information is needed for STK )
         */
        moris_id get_owner() const override;

        // ----------------------------------------------------------------------------

        /**
         * explicitly sets the owner of the edge
         */
        void set_owner( const moris_id aOwner ) override;

        // ----------------------------------------------------------------------------

        /**
         * fills a Vector with pointers to connected vertices
         */
        Vector< mtk::Vertex* > get_vertex_pointers() const override;

        // ----------------------------------------------------------------------------

        void remove_vertex_pointer( moris_index aIndex ) override;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with IDs of connected vertices
         */
        Matrix< IdMat > get_vertex_ids() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with indices of connected vertices
         */
        Matrix< IndexMat > get_vertex_inds() const override;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat > get_vertex_coords() const override;

        // ----------------------------------------------------------------------------

        /**
         * an edge is always a line
         */
        mtk::Geometry_Type get_geometry_type() const override;

        // ----------------------------------------------------------------------------

        mtk::Interpolation_Order get_interpolation_order() const override = 0;

        // ----------------------------------------------------------------------------

        void set_index( const moris_index aIndex ) override;

        // ----------------------------------------------------------------------------

        void set_id( const moris_id aID ) override;

        // ----------------------------------------------------------------------------

        uint get_number_of_elements() const;

        // ----------------------------------------------------------------------------

        Element* get_element( uint aIndex );

        // ----------------------------------------------------------------------------

        uint get_index_on_element( uint aIndex ) const;

        // ----------------------------------------------------------------------------

        Element* get_hmr_leader();

        // ----------------------------------------------------------------------------

        uint get_index_on_leader() const;

        // ----------------------------------------------------------------------------

        virtual const Basis* get_basis( const uint aIndex ) const = 0;

        // ----------------------------------------------------------------------------

        virtual Basis* get_basis( const uint aIndex ) = 0;

        // ----------------------------------------------------------------------------

        bool is_active() const;

        // ----------------------------------------------------------------------------

      private:
        // ----------------------------------------------------------------------------

        void find_leader(
                Mesh_Base*       aMesh,
                Background_Edge* aBackgroundEdge );

        // ----------------------------------------------------------------------------

        virtual void copy_vertex_pointers() = 0;

        // ----------------------------------------------------------------------------

    };    // class Edge

    // ----------------------------------------------------------------------------

}    // namespace moris::hmr
