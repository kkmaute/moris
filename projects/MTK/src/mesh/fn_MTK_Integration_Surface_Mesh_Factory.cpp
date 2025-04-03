/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Surface_Mesh_Data.cpp
 *
 */

#include "fn_MTK_Integration_Surface_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "fn_trans.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"

namespace moris::mtk
{
    Integration_Surface_Mesh_Data::Integration_Surface_Mesh_Data(
            Integration_Mesh const      *aIGMesh,
            const Vector< std::string > &aSideSetNames )
            : Integration_Surface_Mesh_Data( dynamic_cast< Integration_Mesh_DataBase_IG const * >( aIGMesh ), aSideSetNames )
    {
    }

    Integration_Surface_Mesh_Data::Integration_Surface_Mesh_Data(
            Integration_Mesh const     *aIGMesh,
            Vector< Side_Set const * > &aSideSets )
            : mIGMesh( dynamic_cast< Integration_Mesh_DataBase_IG const * >( aIGMesh ) )
    {
        this->initialize_from_side_sets( aSideSets );
        this->initialize_vertex_coordinates();
    }

    Integration_Surface_Mesh_Data::Integration_Surface_Mesh_Data(
            Integration_Mesh_DataBase_IG const *aIGMesh,
            const Vector< std::string >        &aSideSetNames )
            : mIGMesh( aIGMesh )
    {
        Vector< Side_Set const * > tSideSets;

        auto tSideSetFromName = [ &aIGMesh ]( std::string const &aSideSetName ) {
            return dynamic_cast< Side_Set * >( aIGMesh->get_set_by_name( aSideSetName ) );
        };

        std::transform( aSideSetNames.begin(), aSideSetNames.end(), std::back_inserter( tSideSets ), tSideSetFromName );

        this->initialize_from_side_sets( tSideSets );
        this->initialize_vertex_coordinates();
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Side_Set const * >
    Integration_Surface_Mesh_Data::obtain_sidesets_from_names( Integration_Mesh_DataBase_IG const *aIGMesh, const Vector< std::string > &aSideSetNames )
    {
        Vector< Side_Set const * > aSideSets;

        auto aSideSetFromName = [ &aIGMesh ]( std::string const &aSideSetName ) -> Side_Set const * {
            return const_cast< Side_Set const * >( dynamic_cast< Side_Set * >( aIGMesh->get_set_by_name( aSideSetName ) ) );
        };

        std::transform( aSideSetNames.begin(), aSideSetNames.end(), std::back_inserter( aSideSets ), aSideSetFromName );

        return aSideSets;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_from_side_sets( Vector< Side_Set const * > const &aSideSets )
    {
        mSideSets = aSideSets;

        // temporary map to store the neighbors of each vertex since we do not necessarily know the index of the
        // vertex that we want to add as a neighbor at the time of creation.
        map< moris_index, Vector< moris_index > > tTmpNeighborMap;

        // loop over all side sets by name
        for ( auto const &tSideSet : aSideSets )
        {
            this->initialize_side_set( tTmpNeighborMap, dynamic_cast< Set const * >( tSideSet ) );
        }

        // in a last step, the neighbors can actually be correctly assigned since all local indices are known
        this->initialize_neighbors( tTmpNeighborMap );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_side_set(
            map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
            Set const                                 *aSideSet )
    {
        // loop over all clusters that the side set consists of

        moris_index const tNumClustersOnSet = aSideSet->get_clusters_on_set().size();
        mClusterToCellIndices.resize( tNumClustersOnSet, Vector< moris_index >() );
        for ( moris_index tClusterIndex = 0; tClusterIndex < tNumClustersOnSet; tClusterIndex++ )
        {
            Cluster const *tCluster = aSideSet->get_clusters_by_index( tClusterIndex );
            this->initialize_cluster( aTmpNeighborMap, tCluster, tClusterIndex );
        }    // end loop over clusters
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_cluster(
            map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
            Cluster const *const                      &aCluster,
            moris_index                                aClusterIndex )
    {
        Vector< const Cell * > tCells    = aCluster->get_primary_cells_in_cluster();
        Matrix< IdMat >        tCellOrds = aCluster->get_cell_side_ordinals();

        // Cell ordinals define, which side of the cell is actually on the side of the cluster.
        // Each cell should have exactly one edge/facet on the side.
        MORIS_ASSERT( tCells.size() == tCellOrds.size( 1 ), "Number of cells and cell ordinals do not match" );

        // loop over all cells to extract the vertex indices of the ordinals that are actually on the side
        for ( uint i = 0; i < tCells.size(); i++ )
        {
            Cell const *tCurrentCell        = tCells( i );
            int const   tCurrentCellOrdinal = tCellOrds( i );

            this->initialize_cell( aTmpNeighborMap, tCurrentCell, tCurrentCellOrdinal, aClusterIndex );
        }    // end loop over cells
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_cell(
            map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
            const Cell                                *aCell,
            int                                        aCellOrdinal,
            moris_index                                aClusterIndex )
    {
        MORIS_ASSERT( mGlobalToLocalCellIndex.count( aCell->get_index() ) == 0, "Cell added twice to surface mesh" );

        auto const tCurrentLocalCellIndex = static_cast< moris_index >( this->mCellToVertexIndices.size() );

        // local index (on the surface mesh, from 0 to n_surfacemesh), global index (in the integration mesh, arbitrary numbers between 0 and n_igmesh)
        mLocalToGlobalCellIndex.push_back( aCell->get_index() );
        mGlobalToLocalCellIndex[ aCell->get_index() ] = tCurrentLocalCellIndex;

        // one cluster per cell but one cluster can have multiple cells
        mCellToClusterIndices.push_back( aClusterIndex );
        mClusterToCellIndices( aClusterIndex ).push_back( tCurrentLocalCellIndex );

        // prepare the cell to vertex map
        mCellToVertexIndices.push_back( Vector< moris_index >() );

        // side ordinal holds the index of the side of the cell that is actually on the surface
        mCellSideOrdinals.push_back( aCellOrdinal );

        Vector< Vertex const * > tSideVertices = aCell->get_geometric_vertices_on_side_ordinal( aCellOrdinal );

        for ( unsigned int j = 0; j < tSideVertices.size(); j++ )
        {
            Vertex const *tVertex = tSideVertices( j );
            this->initialize_vertex( aTmpNeighborMap, tCurrentLocalCellIndex, tSideVertices, tVertex );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_vertex(
            map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
            moris_index                                aCurrentLocalCellIndex,
            Vector< Vertex const * >                  &aSideVertices,
            Vertex const                              *aVertex )
    {
        moris_index const tVertexIndex             = aVertex->get_index();
        moris_index       tCurrentLocalVertexIndex = 0;
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
            this->mVertexToCellIndices.push_back( Vector< moris_index >() );
        }

        // update the vertex to cell and cell to vertex map for this vertex
        this->mVertexToCellIndices( tCurrentLocalVertexIndex ).push_back( aCurrentLocalCellIndex );
        this->mCellToVertexIndices( aCurrentLocalCellIndex ).push_back( tCurrentLocalVertexIndex );

        for ( auto const &tNeighbor : aSideVertices )
        {    // update neighbors for this vertex for this cell
            if ( tNeighbor != aVertex )
            {
                aTmpNeighborMap[ tVertexIndex ].push_back( tNeighbor->get_index() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_neighbors( map< moris_index, Vector< moris_index > > &aTmpNeighborMap )
    {
        // for each key (vertex in global indices) in the map, the neighbors are assigned.
        mVertexNeighbors.resize( mLocalToGlobalVertexIndex.size() );
        for ( auto const &[ tVertex, tNeighbor ] : aTmpNeighborMap )
        {
            auto const tLocalVertexIndex = mGlobalToLocalVertexIndex[ tVertex ];
            for ( auto const &tNeighborIndex : tNeighbor )
            {
                this->mVertexNeighbors( tLocalVertexIndex ).push_back( mGlobalToLocalVertexIndex[ tNeighborIndex ] );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh_Data::initialize_vertex_coordinates()
    {
        auto const tNumVertices = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        uint const tDim         = mIGMesh->get_spatial_dim();
        mVertexCoordinates.resize( tDim, tNumVertices );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            mVertexCoordinates.set_column( i, mIGMesh->get_node_coordinate( mLocalToGlobalVertexIndex( i ) ) );
        }
    }

}    // namespace moris::mtk
