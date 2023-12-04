/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Surface_Mesh.cpp
 *
 */

#include "cl_MTK_Surface_Mesh.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "fn_trans.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"

namespace moris::mtk
{
    Surface_Mesh::Surface_Mesh( Integration_Mesh *aIGMesh, const moris::Cell< std::string > &aSideSetNames )
            : Surface_Mesh( dynamic_cast< Integration_Mesh_DataBase_IG * >( aIGMesh ), aSideSetNames )
    {
    }


    Surface_Mesh::Surface_Mesh( Integration_Mesh_DataBase_IG *aIGMesh, const moris::Cell< std::string > &aSideSetNames )
            : mIGMesh( aIGMesh )
            , mSideSetNames( aSideSetNames )
    {
        initialize_surface_mesh();
    }

    void Surface_Mesh::initialize_surface_mesh()
    {
        // temporary map to store the neighbors of each vertex since we do not necessarily know the index of the
        // vertex that we want to add as a neighbor at the time of creation.
        map< moris_index, moris::Cell< moris_index > > tTmpNeighborMap;

        // loop over all side sets by name
        for ( auto &tSetName : mSideSetNames )
        {
            Set *tSideSet = mIGMesh->get_set_by_name( tSetName );
            this->initialize_side_set( tTmpNeighborMap, tSideSet );
        }

        // in a last step, the neighbors can actually be correctly assigned since all local indices are known
        this->initialize_neighbors( tTmpNeighborMap );
    }
    void Surface_Mesh::initialize_side_set( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, Set const *aSideSet )
    {
        // loop over all clusters that the side set consists of
        for ( auto &tCluster : aSideSet->get_clusters_on_set() )
        {
            this->initialize_cluster( aTmpNeighborMap, tCluster );
        }    // end loop over clusters
    }

    void Surface_Mesh::initialize_cluster( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, Cluster const *const &aCluster )
    {
        moris::Cell< const Cell * > tCells    = aCluster->get_primary_cells_in_cluster();
        Matrix< IdMat >             tCellOrds = aCluster->get_cell_side_ordinals();

        // Cell ordinals define, which side of the cell is actually on the side of the cluster.
        // Each cell should have exactly one edge/facet on the side.
        MORIS_ASSERT( tCells.size() == tCellOrds.size( 1 ), "Number of cells and cell ordinals do not match" );

        // loop over all cells to extract the vertex indices of the ordinals that are actually on the side
        for ( uint i = 0; i < tCells.size(); i++ )
        {
            Cell const *tCurrentCell        = tCells( i );
            int         tCurrentCellOrdinal = tCellOrds( i );

            this->initialize_cell( aTmpNeighborMap, tCurrentCell, tCurrentCellOrdinal );
        }    // end loop over cells
    }

    void Surface_Mesh::initialize_cell( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, const Cell *aCell, int aCellOrdinal )
    {
        MORIS_ASSERT( mGlobalToLocalCellIndex.count( aCell->get_index() ) == 0, "Cell added twice to surface mesh" );

        auto tCurrentLocalCellIndex = static_cast< moris_index >( this->mCellToVertexIndices.size() );
        mLocalToGlobalCellIndex.push_back( aCell->get_index() );
        mGlobalToLocalCellIndex[ aCell->get_index() ] = tCurrentLocalCellIndex;
        mCellToVertexIndices.push_back( moris::Cell< moris_index >() );
        mCellSideOrdinals.push_back( aCellOrdinal );

        moris::Cell< Vertex const * > tSideVertices = aCell->get_geometric_vertices_on_side_ordinal( aCellOrdinal );

        for ( unsigned int j = 0; j < tSideVertices.size(); j++ )
        {
            Vertex const *tVertex = tSideVertices( j );
            this->initialize_vertex( aTmpNeighborMap, tCurrentLocalCellIndex, tSideVertices, tVertex );
        }
    }

    void Surface_Mesh::initialize_vertex( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, moris_index aCurrentLocalCellIndex, moris::Cell< Vertex const * > &aSideVertices, Vertex const *aVertex )
    {
        moris_index tVertexIndex = aVertex->get_index();
        moris_index tCurrentLocalVertexIndex;
        if ( this->mGlobalToLocalVertexIndex.key_exists( tVertexIndex ) )
        {    // check if the vertex has already been added to the list of vertices. If so, use the local index of the vertex based on the global index.
            tCurrentLocalVertexIndex = this->mGlobalToLocalVertexIndex[ tVertexIndex ];
        }
        else
        {    // if the vertex has not been added to the list of vertices, add it and use the local index based on the current size of the list of vertices
            tCurrentLocalVertexIndex = static_cast< moris_index >( this->mLocalToGlobalVertexIndex.size() );
        }

        if ( this->mGlobalToLocalVertexIndex.count( tVertexIndex ) == 0 )
        {    // check that the vertex has not already been added to the surface mesh
            this->mGlobalToLocalVertexIndex[ tVertexIndex ] = tCurrentLocalVertexIndex;
            this->mLocalToGlobalVertexIndex.push_back( tVertexIndex );
            this->mCellToVertexIndices( aCurrentLocalCellIndex ).push_back( tCurrentLocalVertexIndex );
            this->mVertexToCellIndices.push_back( moris::Cell< moris_index >() );
        }

        // update the vertex to cell map for this vertex
        this->mVertexToCellIndices( tCurrentLocalVertexIndex ).push_back( aCurrentLocalCellIndex );

        for ( auto &tNeighbor : aSideVertices )
        {    // update neighbors for this vertex for this cell
            if ( tNeighbor != aVertex )
            {
                aTmpNeighborMap[ tVertexIndex ].push_back( tNeighbor->get_index() );
            }
        }
    }

    void Surface_Mesh::initialize_neighbors( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap )
    {
        // for each key (vertex in global indices) in the map, the neighbors are assigned.
        mVertexNeighbors.resize( mLocalToGlobalVertexIndex.size() );
        for ( auto &tNeighborPair : aTmpNeighborMap )
        {
            auto tLocalVertexIndex = mGlobalToLocalVertexIndex[ tNeighborPair.first ];
            for ( auto &tNeighborIndex : tNeighborPair.second )
            {
                this->mVertexNeighbors( tLocalVertexIndex ).push_back( mGlobalToLocalVertexIndex[ tNeighborIndex ] );
            }
        }
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates() const
    {
        auto             tNumVertices = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        Matrix< DDRMat > tCoordinates{ tNumVertices, mIGMesh->get_spatial_dim() };
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCoordinates.set_row( i, mIGMesh->get_node_coordinate( mLocalToGlobalVertexIndex( i ) ) );
        }
        if ( mDisplacements.n_rows() > 0 )
        {
            tCoordinates += mDisplacements;
        }
        return tCoordinates;
    }

    moris::Cell< moris::Cell< moris_index > > Surface_Mesh::get_vertex_neighbors() const
    {
        return mVertexNeighbors;
    }

    Matrix< DDRMat > Surface_Mesh::get_facet_normals()
    {
        // only recalculate facet normals if they have not been calculated before
        if ( mFacetNormals.n_rows() == 0 )
        {
            mFacetNormals.resize( mLocalToGlobalCellIndex.size(), mIGMesh->get_spatial_dim() );

            for ( moris::size_t i = 0; i < mLocalToGlobalCellIndex.size(); i++ )
            {
                auto tGlobalCellIndex = mLocalToGlobalCellIndex( i );
                auto tCell            = dynamic_cast< mtk::Cell_DataBase            &>( mIGMesh->get_mtk_cell( tGlobalCellIndex ) );
                auto tSideOrdinal     = mCellSideOrdinals( i );
                auto tNormal          = tCell.compute_outward_side_normal( tSideOrdinal );

                // FIXME: Due to a bug in the compute_outward_side_normal function, the normal, the sign of first component of the normal has to be flipped Until this is not fixed, the following line is used to flip the sign of the first component of the normal
                tNormal( 0, 0 ) *= -1.0;

                mFacetNormals.set_row( i, trans( tNormal ) );
            }
        }

        return mFacetNormals;
    }

    Matrix< DDRMat > Surface_Mesh::get_facet_measure()
    {
        // only recalculate facet measures if they have not been calculated before
        if ( mFacetMeasure.n_rows() == 0 )
        {
            mFacetMeasure.resize( mLocalToGlobalCellIndex.size(), 1 );
            for ( moris::size_t i = 0; i < mLocalToGlobalCellIndex.size(); i++ )
            {
                auto tGlobalCellIndex = mLocalToGlobalCellIndex( i );
                auto tCell            = dynamic_cast< mtk::Cell_DataBase            &>( mIGMesh->get_mtk_cell( tGlobalCellIndex ) );
                auto tSideOrdinal     = mCellSideOrdinals( i );
                auto tMeasure         = tCell.compute_cell_side_measure( tSideOrdinal );
                mFacetMeasure( i )    = tMeasure;
            }
        }

        return mFacetMeasure;
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_normals()
    {
        // only recalculate vertex normals if they have not been calculated before
        if ( mVertexNormals.n_rows() == 0 )
        {
            auto tFacetNormals = this->get_facet_normals();
            auto tFacetMeasure = this->get_facet_measure();
            mVertexNormals.resize( mLocalToGlobalVertexIndex.size(), mIGMesh->get_spatial_dim() );

            auto   tNormal = Matrix< DDRMat >( 1, mIGMesh->get_spatial_dim() );
            size_t tNumNeighbors;
            for ( moris::size_t i = 0; i < mLocalToGlobalVertexIndex.size(); i++ )
            {
                moris::Cell< moris_index > tVertexCellNeighbors = mVertexToCellIndices( i );
                tNumNeighbors                                   = static_cast< moris::size_t >( tVertexCellNeighbors.size() );
                tNormal.fill( 0.0 );    // reset the current normal to zero for each vertex normal calculation

                // compute the normal as the weighted average of the facet normals of the neighboring cells
                for ( moris::size_t j = 0; j < tNumNeighbors; j++ )
                {
                    auto tCellIndex = tVertexCellNeighbors( j );
                    tNormal += tFacetNormals.get_row( tCellIndex ) * tFacetMeasure( tCellIndex );
                }
                mVertexNormals.set_row( i, tNormal / norm( tNormal ) );
            }
        }
        return mVertexNormals;
    }
    moris_index Surface_Mesh::get_global_vertex_index( moris_index aLocalVertexIndex )
    {
        return mLocalToGlobalVertexIndex( aLocalVertexIndex );
    }
    moris_index Surface_Mesh::get_global_cell_index( moris_index aLocalCellIndex )
    {
        return mLocalToGlobalCellIndex( aLocalCellIndex );
    }
    moris_index Surface_Mesh::get_local_vertex_index( moris_index aGlobalVertexIndex )
    {
        return mGlobalToLocalVertexIndex[ aGlobalVertexIndex ];
    }
    moris_index Surface_Mesh::get_local_cell_index( moris_index aGlobalCellIndex )
    {
        return mGlobalToLocalCellIndex[ aGlobalCellIndex ];
    }
    void Surface_Mesh::set_displacement( Matrix< DDRMat > aDisplacements )
    {
        auto nVertices = mLocalToGlobalVertexIndex.size();
        auto nDim      = mIGMesh->get_spatial_dim();
        MORIS_ASSERT( aDisplacements.n_rows() == nVertices, "Number of vertices in displacement matrix does not match number of vertices in mesh" );
        MORIS_ASSERT( aDisplacements.n_cols() == nDim, "Number of dimensions in displacement matrix does not match number of dimensions in mesh" );
        mDisplacements = aDisplacements;
    }


}    // namespace moris::mtk