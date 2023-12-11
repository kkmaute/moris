/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_HMR_Edge.hpp
 *
 */

#ifndef PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_
#define PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_

#include "cl_HMR_Element.hpp"
#include "cl_HMR_Mesh_Base.hpp"
#include "typedefs.hpp"
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
        moris::Cell< Element* > mElements;

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
        virtual ~Edge(){};

        // ----------------------------------------------------------------------------

        /**
         * returns the domain wide id of the edge
         *
         * @return moris_id ID
         */
        moris_id get_id() const;

        // ----------------------------------------------------------------------------

        /**
         * returns the local index of the edge
         *
         * @return moris_index ID
         */
        moris_index get_index() const;

        // ----------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this edge
         */
        virtual uint get_number_of_vertices() const = 0;

        // ----------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this edge
         * ( this information is needed for STK )
         */
        moris_id get_owner() const;

        // ----------------------------------------------------------------------------

        /**
         * explicitly sets the owner of the edge
         */
        void set_owner( const moris_id& aOwner );

        // ----------------------------------------------------------------------------

        /**
         * fills a moris::cell with pointers to connected vertices
         */
        moris::Cell< mtk::Vertex* > get_vertex_pointers() const;

        // ----------------------------------------------------------------------------

        void remove_vertex_pointer( moris_index aIndex );

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with IDs of connected vertices
         */
        Matrix< IdMat > get_vertex_ids() const;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat with indices of connected vertices
         */
        virtual Matrix< IndexMat > get_vertex_inds() const;

        // ----------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat > get_vertex_coords() const;

        // ----------------------------------------------------------------------------

        /**
         * an edge is always a line
         */
        mtk::Geometry_Type get_geometry_type() const;

        // ----------------------------------------------------------------------------

        virtual mtk::Interpolation_Order get_interpolation_order() const = 0;

        // ----------------------------------------------------------------------------

        void set_index( const moris_index& aIndex );

        // ----------------------------------------------------------------------------

        void set_id( const moris_id& aID );

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

        virtual const Basis_Function* get_basis_function( const uint aIndex ) const = 0;

        // ----------------------------------------------------------------------------

        virtual Basis_Function* get_basis_function( const uint aIndex ) = 0;

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

    }; // class Edge

    // ----------------------------------------------------------------------------

} /* namespace moris */

#endif /* PROJECTS_HMR_SRC_CL_HMR_EDGE_HPP_ */
