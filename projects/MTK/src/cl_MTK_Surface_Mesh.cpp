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
            std::cout << "Side set name: " + tSideSet->get_set_name() << std::endl;

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

        //        std::cout << "Cell index: " << aCurrentCell->get_index() << ", ID: " << aCurrentCell->get_id() << " Ord: " << aCurrentCellOrdinal << std::endl;
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

            // check that the vertex has not already been added to the surface mesh
            if ( mGlobalToLocalVertexIndex.count( tVertexIndex ) == 0 )
            {
                auto tCurrentLocalVertexIndex             = static_cast< moris_index >( mLocalToGlobalVertexIndex.size() );
                mGlobalToLocalVertexIndex[ tVertexIndex ] = tCurrentLocalVertexIndex;
                mLocalToGlobalVertexIndex.push_back( tVertexIndex );
                mCellToVertexIndices( aLocalCellIndex ).push_back( tCurrentLocalVertexIndex );

                if ( static_cast< moris_index >( mVertexToCellIndices.size() ) <= tCurrentLocalVertexIndex )
                {
                    mVertexToCellIndices.push_back( moris::Cell< moris_index >() );
                }
                mVertexToCellIndices( tCurrentLocalVertexIndex ).push_back( aLocalCellIndex );

                std::cout << "Vertex: " << tVertex->get_id() << "(" << tVertexIndex << "/" << tCurrentLocalVertexIndex << ")" << std::endl;
            }

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

    Matrix< DDRMat > Surface_Mesh::get_facet_normals() const
    {
        Matrix< DDRMat > tNormals{ mLocalToGlobalCellIndex.size(), mIGMesh->get_spatial_dim() };
        for ( moris::size_t i = 0; i < mLocalToGlobalCellIndex.size(); i++ )
        {
            auto tGlobalCellIndex = mLocalToGlobalCellIndex( i );
            auto tCell            = dynamic_cast< mtk::Cell_DataBase            &>( mIGMesh->get_mtk_cell( tGlobalCellIndex ) );
            //            auto tSideOrdinal     = mCellSideOrdinals( i );
            for ( int j = 0; j < 3; j++ )
            {
                auto tNormal   = tCell.compute_outward_side_normal( j );
                auto tVertices = tCell.get_vertices_on_side_ordinal( j );

                std::cout << "Cell: " << tCell.get_id() << std::endl;
                std::cout << "Ordinal " << j << " between vertices "
                          << tVertices( 0 )->get_id() << " and "
                          << tVertices( 1 )->get_id() << std::endl;
                std::cout << "Coordinates ("
                          << tVertices( 0 )->get_coords()( 0 ) << ", "
                          << tVertices( 0 )->get_coords()( 1 ) << ") and ("
                          << tVertices( 1 )->get_coords()( 0 ) << ", "
                          << tVertices( 1 )->get_coords()( 1 ) << ")" << std::endl;
                std::cout << "Normal: "
                          << tNormal( 0 ) << " "
                          << tNormal( 1 ) << std::endl;
                std::cout << " -------------------------------- \n"
                          << std::endl;
            }

            //            tNormals.set_row( i, trans( tNormal ) );
        }

        return tNormals;
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_normals() const
    {
        return Matrix< DDRMat >();
    }
}    // namespace moris::mtk