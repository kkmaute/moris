/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell_Proxy.hpp
 *
 */

#ifndef PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_
#define PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Cell_Cluster.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    class Cell_Proxy : public Cell
    {
        //------------------------------------------------------------------------------

      public:
        moris::moris_id    mId    = gNoID;
        moris::moris_index mIndex = gNoIndex;
        moris::moris_id    mOwner = gNoID;

        Vector< Vertex * >       mVertices;
        moris::uint              mSpatialDim   = 3;
        enum Geometry_Type       mGeometryType = Geometry_Type::UNDEFINED;
        enum Interpolation_Order mInterpOrder  = Interpolation_Order::UNDEFINED;
        enum Integration_Order   mIntegOrder   = Integration_Order::UNDEFINED;
        moris::mtk::Cell_Info   *mCellInfo     = nullptr;

        //------------------------------------------------------------------------------

        Cell_Proxy() {};

        //------------------------------------------------------------------------------

        ~Cell_Proxy() override {};

        //------------------------------------------------------------------------------

        /**
         * returns the domain wide id of the cell
         *
         * @return moris_id ID
         */
        moris_id
        get_id() const override
        {
            MORIS_ASSERT( mId != gNoID, "Cell ID not initialized" );

            return mId;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the local index of the cell
         *
         * @return moris_index ID
         */
        moris_index
        get_index() const override
        {
            MORIS_ASSERT( mIndex != gNoIndex, "Cell index not initialized" );

            return mIndex;
        }

        //------------------------------------------------------------------------------

        /**
         * tells how many vertices are connected to this cell
         */
        uint
        get_number_of_vertices() const override
        {
            return mVertices.size();
        }

        //------------------------------------------------------------------------------

        /**
         * returns the proc id of the owner of this cell
         * ( this information is needed for STK )
         */
        moris_id
        get_owner() const override
        {
            MORIS_ASSERT( mOwner != gNoID, "Cell ownership not initialized" );

            return mOwner;
        }

        //------------------------------------------------------------------------------

        /**
         * fills a Vector with pointers to connected vertices
         */
        Vector< Vertex * >
        get_vertex_pointers() const override
        {
            return mVertices;
        }

        //------------------------------------------------------------------------------

        // TODO MESHCLEANUP
        void
        remove_vertex_pointer( moris_index aIndex ) override
        {
            // do nothing for now
            // std::cout<<"In MTK Cell Proxy"<<std::endl;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with IDs of connected vertices
         */
        Matrix< IdMat >
        get_vertex_ids() const override
        {
            Matrix< IdMat > tVertexIds( mVertices.size(), 1 );

            for ( moris::uint i = 0; i < this->get_number_of_vertices(); i++ )
            {
                tVertexIds( i ) = mVertices( i )->get_id();
            }

            return tVertexIds;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat with indices of connected vertices
         */
        Matrix< IndexMat >
        get_vertex_inds() const override
        {
            Matrix< IndexMat > tVertexInd( mVertices.size(), 1 );

            for ( moris::uint i = 0; i < this->get_number_of_vertices(); i++ )
            {
                tVertexInd( i ) = mVertices( i )->get_index();
            }

            return tVertexInd;
        }

        //------------------------------------------------------------------------------

        /**
         * returns a Mat of dimension
         * < number of vertices * number of dimensions >
         */
        Matrix< DDRMat >
        get_vertex_coords() const override
        {
            size_t           tNumVertices = this->get_number_of_vertices();
            Matrix< DDRMat > tVertexCoords( tNumVertices, mSpatialDim );
            for ( size_t i = 0; i < tNumVertices; i++ )
            {
                tVertexCoords.set_row( i, mVertices( i )->get_coords() );
            }
            return tVertexCoords;
        }

        //------------------------------------------------------------------------------

        /**
         * returns an enum that defines the geometry type of the element
         */
        Geometry_Type
        get_geometry_type() const override
        {
            return mGeometryType;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the order of the element
         */
        Interpolation_Order
        get_interpolation_order() const override
        {
            return mInterpOrder;
        }

        //------------------------------------------------------------------------------

        /**
         * returns the integration order of the element
         */
        Integration_Order
        get_integration_order() const override
        {
            return mIntegOrder;
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cell_measure() const override
        {
            return mCellInfo->compute_cell_size( this );
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cell_measure_deriv( uint aLocalVertexID, uint aDirection ) const override
        {
            return mCellInfo->compute_cell_size_deriv_general( this, aLocalVertexID, aDirection );
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cell_side_measure( moris_index const &aSideOrdinal ) const override
        {
            return mCellInfo->compute_cell_side_size( this, aSideOrdinal );
        }

        //------------------------------------------------------------------------------

        moris::real
        compute_cell_side_measure_deriv(
                moris_index const &aCellSideOrd,
                uint               aLocalVertexID,
                uint               aDirection ) const override
        {
            MORIS_ERROR( 0, "compute_cell_side_measure_deriv - Not implemented." );
            return 0.0;
        }

        //------------------------------------------------------------------------------

    };    // class Cell_Proxy

    //------------------------------------------------------------------------------

}    // namespace moris::mtk

#endif /* PROJECTS_MTK_TEST_CL_MTK_CELL_PROXY_HPP_ */
