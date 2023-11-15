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
    {
        initialize_surface_mesh( aSideSetNames );
    }

    void Surface_Mesh::initialize_surface_mesh( moris::Cell< std::string > const &aSideSetNames )
    {
        // loop over all side sets by name
        for ( auto &tSetName : aSideSetNames )
        {
            Set *tSideSet = mIGMesh->get_set_by_name( tSetName );

            // temporary map to store the neighbors of each vertex since we do not necessarily know the index of the
            // vertex that we want to add as a neighbor at the time of creation.
            moris::map< moris_index, moris::Cell< moris_index > > tTmpNeighborMap;

            // loop over all clusters that the side set consists of
            for ( auto &tCluster : tSideSet->get_clusters_on_set() )
            {
                moris::Cell< const moris::mtk::Cell * > tCells    = tCluster->get_primary_cells_in_cluster();
                moris::Matrix< moris::IdMat >           tCellOrds = tCluster->get_cell_side_ordinals();

                // Cell ordinals define, which side of the cell is actually on the side of the cluster.
                // Each cell should have exactly one edge/facet on the side.
                MORIS_ASSERT( tCells.size() == tCellOrds.size( 1 ), "Number of cells and cell ordinals do not match" );

                // loop over all cells to extract the vertex indices of the ordinals that are actually on the side
                for ( moris::uint i = 0; i < tCells.size(); i++ )
                {
                    Cell const *tCurrentCell           = tCells( i );
                    int         tCurrentCellOrdinal    = tCellOrds( i );
                    auto        tCurrentLocalCellIndex = static_cast< moris_index >( mCellToVertexIndices.size() );

                    initialize_cell( tCurrentCell, tCurrentCellOrdinal );
                    initialize_vertices( tTmpNeighborMap, tCurrentCell, tCurrentLocalCellIndex, tCurrentCellOrdinal );
                }    // end loop over cells
            }        // end loop over clusters

            // in a last step, the neighbors can actually be correctly assigned since all local indices are known
            this->initialize_neighbors( tTmpNeighborMap );
        }    // end loop over side sets
    }

    void Surface_Mesh::initialize_cell( const Cell *aCurrentCell, int aCurrentCellOrdinal )
    {
        MORIS_ASSERT( mGlobalToLocalCellIndex.count( aCurrentCell->get_index() ) == 0, "Cell added twice to surface mesh" );

        mLocalToGlobalCellIndex.push_back( aCurrentCell->get_index() );
        mGlobalToLocalCellIndex[ aCurrentCell->get_index() ] = static_cast< moris_index >( mLocalToGlobalCellIndex.size() ) - 1;
        mCellToVertexIndices.push_back( moris::Cell< moris_index >() );
        mCellSideOrdinals.push_back( aCurrentCellOrdinal );
    }

    void Surface_Mesh::initialize_vertices( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, const Cell *aCell, moris_index aLocalCellIndex, int aCellOrdinal )
    {
        moris::Cell< Vertex const * > tSideVertices = aCell->get_geometric_vertices_on_side_ordinal( aCellOrdinal );

        for ( unsigned int j = 0; j < tSideVertices.size(); j++ )
        {
            Vertex const *tVertex      = tSideVertices( j );
            moris_index   tVertexIndex = tVertex->get_index();

            moris_index tCurrentLocalVertexIndex;

            // check if the vertex has already been added to the list of vertices. If so, use the local index of the vertex based on the global index.
            if ( mGlobalToLocalVertexIndex.key_exists( tVertexIndex ) )
            {
                tCurrentLocalVertexIndex = mGlobalToLocalVertexIndex[ tVertexIndex ];
            }
            else
            {
                // if the vertex has not been added to the list of vertices, add it and use the local index based on the current size of the list of vertices
                tCurrentLocalVertexIndex = static_cast< moris_index >( mLocalToGlobalVertexIndex.size() );
            }

            // check that the vertex has not already been added to the surface mesh
            if ( mGlobalToLocalVertexIndex.count( tVertexIndex ) == 0 )
            {
                mGlobalToLocalVertexIndex[ tVertexIndex ] = tCurrentLocalVertexIndex;
                mLocalToGlobalVertexIndex.push_back( tVertexIndex );
                mCellToVertexIndices( aLocalCellIndex ).push_back( tCurrentLocalVertexIndex );
                mVertexToCellIndices.push_back( moris::Cell< moris_index >() );
            }

            // update the vertex to cell map for this vertex
            mVertexToCellIndices( tCurrentLocalVertexIndex ).push_back( aLocalCellIndex );

            // update neighbors for this vertex for this cell
            for ( auto &tNeighbor : tSideVertices )
            {
                if ( tNeighbor != tVertex )
                {
                    aTmpNeighborMap[ tVertexIndex ].push_back( tNeighbor->get_index() );
                }
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


}    // namespace moris::mtk