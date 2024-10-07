/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Cell.cpp
 *
 */

#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Cell_Info.hpp"
#include "cl_Matrix.hpp"
#include "fn_trans.hpp"
#include "fn_cross.hpp"
#include "fn_norm.hpp"
#include <typeinfo>
#include <utility>

namespace moris::mtk
{
    //------------------------------------------------------------------------------
    Cell::Cell( moris_id                             aCellId,
            moris_index                              aCellIndex,
            moris_id                                 aCellOwner,
            std::shared_ptr< moris::mtk::Cell_Info > aCellInfo )
            : mCellInfo( std::move( aCellInfo ) )
            , mCellId( aCellId )
            , mCellIndex( aCellIndex )
            , mCellOwner( aCellOwner )
    {
    }

    //------------------------------------------------------------------------------
    Cell_Info const *
    Cell::get_cell_info() const
    {
        MORIS_ASSERT( mCellInfo != nullptr, "Cell info not set" );
        return mCellInfo.get();
    }
    //------------------------------------------------------------------------------

    std::shared_ptr< mtk::Cell_Info >
    Cell::get_cell_info_sp() const
    {
        MORIS_ASSERT( mCellInfo != nullptr, "Cell info not set" );
        return mCellInfo;
    }

    //------------------------------------------------------------------------------

    void
    Cell::set_mtk_cell_info( std::shared_ptr< moris::mtk::Cell_Info > aCellInfo )
    {
        mCellInfo = std::move( aCellInfo );
    }

    //------------------------------------------------------------------------------

    moris_id
    Cell::get_id() const
    {
        return mCellId;
    }

    void
    Cell::set_id( const moris_id aId )
    {
        mCellId = aId;
    }

    //------------------------------------------------------------------------------

    void
    Cell::set_index( const moris_index aIndex )
    {
        mCellIndex = aIndex;
    }

    //------------------------------------------------------------------------------

    void
    Cell::set_owner( const moris_id aOwner )
    {
        mCellOwner = aOwner;
    }

    //------------------------------------------------------------------------------

    moris_index
    Cell::get_index() const
    {
        return mCellIndex;
    }

    //------------------------------------------------------------------------------

    uint
    Cell::get_level() const
    {
        return 0;
    }

    //------------------------------------------------------------------------------

    moris_id
    Cell::get_owner() const
    {
        return mCellOwner;
    }

    //------------------------------------------------------------------------------

    uint
    Cell::get_number_of_vertices() const
    {
        return this->get_cell_info()->get_num_verts();
    }
    //------------------------------------------------------------------------------

    uint
    Cell::get_number_of_facets() const
    {
        return this->get_cell_info()->get_num_facets();
    }

    //------------------------------------------------------------------------------

    uint
    Cell::get_number_of_edges() const
    {
        return this->get_cell_info()->get_num_edges();
    }
    //------------------------------------------------------------------------------

    Matrix< IdMat >
    Cell::get_vertex_ids() const
    {
        uint              tNumVertices = this->get_number_of_vertices();
        Vector< Vertex* > tVertices    = this->get_vertex_pointers();

        Matrix< IdMat > tVertexIds( 1, tNumVertices );
        for ( uint i = 0; i < tNumVertices; i++ )
        {
            tVertexIds( i ) = tVertices( i )->get_id();
        }
        return tVertexIds;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cell::get_vertex_inds() const
    {
        uint              tNumVertices = this->get_number_of_vertices();
        Vector< Vertex* > tVertices    = this->get_vertex_pointers();

        Matrix< IdMat > tVertexInds( 1, tNumVertices );
        for ( uint i = 0; i < tNumVertices; i++ )
        {
            tVertexInds( i ) = tVertices( i )->get_index();
        }
        return tVertexInds;
    }

    //------------------------------------------------------------------------------

    void
    Cell::remove_vertex( moris_index aIndex )
    {
        this->remove_vertex_pointer( aIndex );
    }

    //------------------------------------------------------------------------------

    bool
    Cell::check_unique_vertex_inds() const
    {
        uint                      tNumVertices = this->get_number_of_vertices();
        moris::Matrix< IndexMat > tVertexInds  = this->get_vertex_inds();
        for ( uint i = 0; i < tNumVertices; i++ )
        {
            for ( uint j = 0; j < i; j++ )
            {
                if ( tVertexInds( i ) == tVertexInds( j ) )
                {
                    return false;
                }
            }
        }
        return true;
    }

    //------------------------------------------------------------------------------

    Matrix< IndexMat >
    Cell::get_vertex_owners() const
    {
        uint              tNumVertices = this->get_number_of_vertices();
        Vector< Vertex* > tVertices    = this->get_vertex_pointers();

        Matrix< IdMat > tVertexOwners( 1, tNumVertices );
        for ( uint i = 0; i < tNumVertices; i++ )
        {
            tVertexOwners( i ) = tVertices( i )->get_owner();
        }
        return tVertexOwners;
    }

    //------------------------------------------------------------------------------

    Vector< mtk::Vertex_Interpolation* >
    Cell::get_vertex_interpolations( const uint aOrder ) const
    {
        uint                                 tNumVerts       = this->get_number_of_vertices();
        Vector< mtk::Vertex* >               tVertexPointers = this->get_vertex_pointers();
        Vector< mtk::Vertex_Interpolation* > tVertexInterp( tNumVerts );

        for ( moris::uint i = 0; i < tNumVerts; i++ )
        {
            tVertexInterp( i ) = tVertexPointers( i )->get_interpolation( aOrder );
        }

        return tVertexInterp;
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Cell::get_vertices_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
        Vector< moris::mtk::Vertex* > tVertices = this->get_vertex_pointers();

        moris::Matrix< moris::IndexMat > tNodeOrdsOnSide = this->get_cell_info()->get_node_to_facet_map( aSideOrdinal );

        Vector< moris::mtk::Vertex const * > tVerticesOnSide( tNodeOrdsOnSide.numel() );
        for ( moris::uint i = 0; i < tNodeOrdsOnSide.numel(); i++ )
        {
            tVerticesOnSide( i ) = tVertices( tNodeOrdsOnSide( i ) );
        }
        return tVerticesOnSide;
    }

    //------------------------------------------------------------------------------

    Vector< moris::mtk::Vertex const * >
    Cell::get_geometric_vertices_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
        MORIS_ASSERT( mCellInfo != nullptr, "Cell info null ptr" );

        Vector< moris::mtk::Vertex* > tVertices = this->get_vertex_pointers();

        // get vertex ordinals
        moris::Matrix< moris::IndexMat > tGeometricVertOrdsOnFacet = this->get_cell_info()->get_geometric_node_to_facet_map( aSideOrdinal );

        // allocate cell of vertices
        Vector< moris::mtk::Vertex const * > tVerticesOnSide( tGeometricVertOrdsOnFacet.numel() );

        for ( moris::uint i = 0; i < tGeometricVertOrdsOnFacet.numel(); i++ )
        {
            tVerticesOnSide( i ) = tVertices( tGeometricVertOrdsOnFacet( i ) );
        }

        return tVerticesOnSide;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell::get_cell_physical_coords_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
        // FIXME: Add assert to check side ordinal

        // get the vertex pointers on the side
        Vector< moris::mtk::Vertex const * > tVerticesOnSide = this->get_vertices_on_side_ordinal( aSideOrdinal );

        // allocate output coords (note we do not know the spatial dimension at this time)
        moris::Matrix< moris::DDRMat > tVertexPhysCoords( 0, 0 );

        // iterate through vertices and collect local coordinates
        for ( moris::uint i = 0; i < tVerticesOnSide.size(); i++ )
        {
            moris::Matrix< moris::DDRMat > tVertexCoord = tVerticesOnSide( i )->get_coords();

            if ( i == 0 )
            {
                MORIS_ASSERT( isrow( tVertexCoord ), "Default implementation assumes row based coordinates" );
                tVertexPhysCoords.resize( tVerticesOnSide.size(), tVertexCoord.numel() );
            }

            tVertexPhysCoords.get_row( i ) = tVertexCoord.get_row( 0 );
        }

        return tVertexPhysCoords;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< IndexMat >
    Cell::get_vertices_ind_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
        Vector< moris::mtk::Vertex const * > tVertices = this->get_vertices_on_side_ordinal( aSideOrdinal );

        uint tNumVertices = tVertices.size();

        Matrix< IndexMat > tVertexInd( 1, tNumVertices );

        for ( uint i = 0; i < tNumVertices; i++ )
        {
            tVertexInd( 0, i ) = tVertices( i )->get_index();
        }
        return tVertexInd;
    }

    //------------------------------------------------------------------------------

    moris_index
    Cell::get_vertex_ordinal_wrt_cell( moris_index const & aVertexIndex ) const
    {
        // get all vertex indices on current mtk::cell
        Matrix< IndexMat > tVertexInds = this->get_vertex_inds();

        // loop through all indices on cell and check if vertex index lives on this cell
        for ( moris::uint i = 0; i < tVertexInds.numel(); i++ )
        {
            if ( tVertexInds( i ) == aVertexIndex )
            {
                return (moris_index)i;
            }
        }

        // throw error if vertex is not found on mtk::cell
        MORIS_ERROR( false, "mtk::Cell::get_vertex_ordinal_wrt_cell() - Vertex not attached to cell." );
        return MORIS_INDEX_MAX;
    }

    //------------------------------------------------------------------------------

    Geometry_Type
    Cell::get_geometry_type() const
    {
        return this->get_cell_info()->get_cell_geometry();
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell::compute_outward_side_normal( moris::moris_index aSideOrdinal ) const
    {
        // get the vertex coordinates
        moris::Matrix< moris::DDRMat > tVertexCoords = this->get_vertex_coords();

        // Get the nodes which need to be used to compute normal
        moris::Matrix< moris::IndexMat > tEdgeNodesForNormal = this->get_cell_info()->get_node_map_outward_normal( aSideOrdinal );

        // initialize the output normal
        Matrix< DDRMat > tOutwardNormal;

        // for 2D just get the normal for an edge
        if ( tVertexCoords.n_cols() == 2 )
        {
            moris_index tFirstNode  = tEdgeNodesForNormal( 0 );
            moris_index tSecondNode = tEdgeNodesForNormal( 1 );
            real        tX1         = tVertexCoords( tFirstNode, 0 );
            real        tX2         = tVertexCoords( tSecondNode, 0 );
            real        tY1         = tVertexCoords( tFirstNode, 1 );
            real        tY2         = tVertexCoords( tSecondNode, 1 );
            tOutwardNormal          = { { tY2 - tY1 }, { tX2 - tX1 } };
        }

        // for 3D get the normal to the plane spanned by two edges
        else
        {
            // Get vector along these edges
            moris::Matrix< moris::DDRMat > tEdge0Vector( tVertexCoords.numel(), 1 );
            moris::Matrix< moris::DDRMat > tEdge1Vector( tVertexCoords.numel(), 1 );

            // Get vector along these edges
            tEdge0Vector = moris::linalg_internal::trans( tVertexCoords.get_row( tEdgeNodesForNormal( 1, 0 ) ) - tVertexCoords.get_row( tEdgeNodesForNormal( 0, 0 ) ) );
            tEdge1Vector = moris::linalg_internal::trans( tVertexCoords.get_row( tEdgeNodesForNormal( 1, 1 ) ) - tVertexCoords.get_row( tEdgeNodesForNormal( 0, 1 ) ) );

            // Take the cross product to get the normal
            tOutwardNormal = moris::cross( tEdge0Vector, tEdge1Vector );
        }

        // Normalize
        Matrix< DDRMat > tUnitOutwardNormal = tOutwardNormal / moris::norm( tOutwardNormal );

        return tUnitOutwardNormal;
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell::compute_cell_measure() const
    {
        return this->get_cell_info()->compute_cell_size( this );
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell::compute_cell_measure_deriv(
            uint aLocalVertexID,
            uint aDirection ) const
    {
        return this->get_cell_info()->compute_cell_size_deriv( this, aLocalVertexID, aDirection );
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell::compute_cell_side_measure( moris_index const & aCellSideOrd ) const
    {
        return this->get_cell_info()->compute_cell_side_size( this, aCellSideOrd );
    }

    //------------------------------------------------------------------------------

    Interpolation_Order
    Cell::get_interpolation_order() const
    {
        return this->get_cell_info()->get_cell_interpolation_order();
    }

    //------------------------------------------------------------------------------

    mtk::Integration_Order
    Cell::get_integration_order() const
    {
        return this->get_cell_info()->get_cell_integration_order();
    }

    //------------------------------------------------------------------------------

    moris::real
    Cell::compute_cell_side_measure_deriv(
            moris_index const & aCellSideOrd,
            uint                aLocalVertexID,
            uint                aDirection ) const
    {
        return this->get_cell_info()->compute_cell_side_size_deriv( this, aCellSideOrd, aLocalVertexID, aDirection );
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell::get_cell_geometric_coords_on_side_ordinal( moris::moris_index aSideOrdinal ) const
    {
        MORIS_ASSERT( mCellInfo != nullptr, "Cell info null ptr" );

        // Get all the vertices of the cell
        Vector< moris::mtk::Vertex* > tVertices = this->get_vertex_pointers();

        // Get vertex ordinals on the facet
        moris::Matrix< moris::IndexMat > tGeometricVertOrdsOnFacet = this->get_cell_info()->get_geometric_node_to_facet_map( aSideOrdinal );

        // Allocate output coords (note we do not know the spatial dimension at this time)
        moris::Matrix< moris::DDRMat > tVertexPhysCoords( 0, 0 );

        for ( moris::uint i = 0; i < tGeometricVertOrdsOnFacet.numel(); i++ )
        {
            moris::Matrix< moris::DDRMat > tVertexCoord = tVertices( tGeometricVertOrdsOnFacet( i ) )->get_coords();

            if ( i == 0 )
            {
                MORIS_ASSERT( isrow( tVertexCoord ), "Default implementation assumes row based coordinates" );
                tVertexPhysCoords.resize( tGeometricVertOrdsOnFacet.numel(), tVertexCoord.numel() );
            }

            tVertexPhysCoords.get_row( i ) = tVertexCoord.get_row( 0 );
        }

        return tVertexPhysCoords;
    }

    //------------------------------------------------------------------------------

    moris::Matrix< moris::DDRMat >
    Cell::compute_cell_centroid() const
    {
        // get all vertices of cell
        Vector< Vertex* > tVertices = this->get_vertex_pointers();

        // get coordinate of first vertex
        Matrix< DDRMat > tCentroid = tVertices( 0 )->get_coords();

        // add coordinates of remaining vertices
        for ( uint iI = 1; iI < tVertices.size(); ++iI )
        {
            tCentroid += tVertices( iI )->get_coords();
        }

        // compute and return coordinate average
        return 1.0 / tVertices.size() * tCentroid;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::mtk
