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
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    Surface_Mesh::Surface_Mesh( Integration_Mesh *aIGMesh, const moris::Cell< std::string > &aSideSetNames )
            : Surface_Mesh( dynamic_cast< Integration_Mesh_DataBase_IG * >( aIGMesh ), aSideSetNames )
    {
    }

    Surface_Mesh::Surface_Mesh( Integration_Mesh *aIGMesh, moris::Cell< Side_Set const * > &aSideSets )
            : mIGMesh( dynamic_cast< Integration_Mesh_DataBase_IG * >( aIGMesh ) )
    {
        this->initialize_from_side_sets( aSideSets );
    }

    Surface_Mesh::Surface_Mesh( Integration_Mesh_DataBase_IG *aIGMesh, const moris::Cell< std::string > &aSideSetNames )
    {
        moris::Cell< Side_Set const * > tSideSets;

        auto tSideSetFromName = [ &aIGMesh ]( std::string const &aSideSetName ) {
            return dynamic_cast< Side_Set * >( aIGMesh->get_set_by_name( aSideSetName ) );
        };
        std::transform( aSideSetNames.begin(), aSideSetNames.end(), std::back_inserter( tSideSets ), tSideSetFromName );

        initialize_from_side_sets( tSideSets );
    }

    void Surface_Mesh::initialize_from_side_sets( moris::Cell< Side_Set const * > const &aSideSets )
    {
        mSideSets = aSideSets;

        // temporary map to store the neighbors of each vertex since we do not necessarily know the index of the
        // vertex that we want to add as a neighbor at the time of creation.
        map< moris_index, moris::Cell< moris_index > > tTmpNeighborMap;

        // loop over all side sets by name
        for ( auto const &tSideSet : aSideSets )
        {
            this->initialize_side_set( tTmpNeighborMap, dynamic_cast< Set const * >( tSideSet ) );
        }

        // in a last step, the neighbors can actually be correctly assigned since all local indices are known
        this->initialize_neighbors( tTmpNeighborMap );
        this->initialize_vertex_coordinates();
        this->initialize_facet_normals();
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    void Surface_Mesh::initialize_side_set( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, Set const *aSideSet )
    {
        // loop over all clusters that the side set consists of

        moris_index const tNumClustersOnSet = aSideSet->get_clusters_on_set().size();
        mClusterToCellIndices.resize( tNumClustersOnSet, moris::Cell< moris_index >() );
        for ( moris_index tClusterIndex = 0; tClusterIndex < tNumClustersOnSet; tClusterIndex++ )
        {
            Cluster const *tCluster = aSideSet->get_clusters_by_index( tClusterIndex );
            this->initialize_cluster( aTmpNeighborMap, tCluster, tClusterIndex );
        }    // end loop over clusters
    }

    void Surface_Mesh::initialize_cluster( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, Cluster const *const &aCluster, moris_index aClusterIndex )
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
            int const   tCurrentCellOrdinal = tCellOrds( i );

            this->initialize_cell( aTmpNeighborMap, tCurrentCell, tCurrentCellOrdinal, aClusterIndex );
        }    // end loop over cells
    }

    void Surface_Mesh::initialize_cell( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap, const Cell *aCell, int aCellOrdinal, moris_index aClusterIndex )
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
        mCellToVertexIndices.push_back( moris::Cell< moris_index >() );

        // side ordinal holds the index of the side of the cell that is actually on the surface
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
            this->mVertexToCellIndices.push_back( moris::Cell< moris_index >() );
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

    void Surface_Mesh::initialize_neighbors( map< moris_index, moris::Cell< moris_index > > &aTmpNeighborMap )
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

    void Surface_Mesh::initialize_vertex_coordinates()
    {
        auto const tNumVertices = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        uint const tDim         = mIGMesh->get_spatial_dim();
        mVertexCoordinates.resize( tDim, tNumVertices );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            mVertexCoordinates.set_column( i, mIGMesh->get_node_coordinate( mLocalToGlobalVertexIndex( i ) ) );
        }
    }

    void Surface_Mesh::initialize_facet_normals()
    {
        auto const tNumCells = static_cast< moris::size_t >( mLocalToGlobalCellIndex.size() );
        uint const tDim      = mIGMesh->get_spatial_dim();
        mFacetNormals.resize( tDim, tNumCells );

        for ( moris::size_t i = 0; i < tNumCells; i++ )
        {
            auto const tGlobalCellIndex = mLocalToGlobalCellIndex( i );
            auto       tCell            = dynamic_cast< mtk::Cell_DataBase                  &>( mIGMesh->get_mtk_cell( tGlobalCellIndex ) );
            auto const tSideOrdinal     = mCellSideOrdinals( i );
            auto       tNormal          = tCell.compute_outward_side_normal( tSideOrdinal );

            // FIXME: Due to a bug in the compute_outward_side_normal function, the normal, the sign of first component of the normal has to be flipped Until this is not fixed, the following line is used to flip the sign of the first component of the normal
            tNormal( 0, 0 ) *= -1.0;

            mFacetNormals.set_column( i, tNormal );
        }
    }

    void Surface_Mesh::initialize_facet_measure()
    {
        auto const tNumCells = static_cast< moris::size_t >( mLocalToGlobalCellIndex.size() );
        mFacetMeasure.resize( tNumCells, 1 );
        for ( moris::size_t i = 0; i < tNumCells; i++ )
        {
            int const   tGlobalCellIndex = mLocalToGlobalCellIndex( i );
            int         tSideOrdinal     = mCellSideOrdinals( i );
            Cell const &tCell            = mIGMesh->get_mtk_cell( tGlobalCellIndex );
            real const  tMeasure         = tCell.compute_cell_side_measure( tSideOrdinal );
            mFacetMeasure( i )           = tMeasure;
        }
    }

    void Surface_Mesh::initialize_vertex_normals()
    {
        auto const                    tNumVertices  = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        uint const                    tDim          = mIGMesh->get_spatial_dim();
        Matrix< arma::Mat< double > > tFacetNormals = this->get_facet_normals();
        Matrix< arma::Mat< double > > tFacetMeasure = this->get_facet_measure();
        mVertexNormals.resize( tDim, tNumVertices );

        auto tNormal = Matrix< DDRMat >( tDim, 1 );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            moris::Cell< moris_index > tVertexCellNeighbors = mVertexToCellIndices( i );
            auto const                 tNumNeighbors        = static_cast< moris::size_t >( tVertexCellNeighbors.size() );
            tNormal.fill( 0.0 );    // reset the current normal to zero for each vertex normal calculation

            // compute the normal as the weighted average of the facet normals of the neighboring cells
            for ( moris::size_t j = 0; j < tNumNeighbors; j++ )
            {
                int const tCellIndex = tVertexCellNeighbors( j );
                tNormal += tFacetNormals.get_column( tCellIndex ) * tFacetMeasure( tCellIndex );
            }
            mVertexNormals.set_column( i, tNormal / norm( tNormal ) );
        }
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates() const
    {
        if ( mDisplacements.n_cols() > 0 )
        {
            return mVertexCoordinates + mDisplacements;
        }
        return mVertexCoordinates;
    }

    moris::Cell< moris::Cell< moris_index > > Surface_Mesh::get_vertex_neighbors() const
    {
        return mVertexNeighbors;
    }

    moris::Cell< moris_index > Surface_Mesh::get_vertex_neighbors( moris_index aLocalVertexIndex ) const
    {
        MORIS_ASSERT( aLocalVertexIndex < static_cast< moris_index >( mVertexNeighbors.size() ), "Vertex index out of bounds" );
        return mVertexNeighbors( aLocalVertexIndex );
    }

    Matrix< DDRMat > Surface_Mesh::get_facet_normals() const
    {
        return mFacetNormals;
    }

    Matrix< DDRMat > Surface_Mesh::get_facet_measure() const
    {
        return mFacetMeasure;
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_normals() const
    {
        return mVertexNormals;
    }

    moris_index Surface_Mesh::get_global_vertex_index( moris_index aLocalVertexIndex ) const
    {
        return mLocalToGlobalVertexIndex( aLocalVertexIndex );
    }

    moris_index Surface_Mesh::get_global_cell_index( moris_index aLocalCellIndex ) const
    {
        return mLocalToGlobalCellIndex( aLocalCellIndex );
    }

    moris_index Surface_Mesh::get_local_vertex_index( moris_index aGlobalVertexIndex ) const
    {
        return mGlobalToLocalVertexIndex.at( aGlobalVertexIndex );
    }

    moris_index Surface_Mesh::get_local_cell_index( moris_index aGlobalCellIndex ) const
    {
        return mGlobalToLocalCellIndex.at( aGlobalCellIndex );
    }

    void Surface_Mesh::set_displacement( Matrix< DDRMat > const &aDisplacements )
    {
        MORIS_ASSERT( aDisplacements.n_rows() == mLocalToGlobalVertexIndex.size(), "Number of vertices in displacement matrix does not match number of vertices in mesh" );
        MORIS_ASSERT( aDisplacements.n_cols() == mIGMesh->get_spatial_dim(), "Number of dimensions in displacement matrix does not match number of dimensions in mesh" );
        mDisplacements = aDisplacements;
    }

    moris_index Surface_Mesh::get_cluster_of_cell( moris_index aLocalCellIndex ) const
    {
        return mCellToClusterIndices( aLocalCellIndex );
    }

    moris::Cell< moris_index > Surface_Mesh::get_vertices_of_cell( moris_index aLocalCellIndex ) const
    {
        return mCellToVertexIndices( aLocalCellIndex );
    }

    moris::Cell< moris_index > Surface_Mesh::get_cells_of_vertex( moris_index aLocalVertexIndex ) const
    {
        return mVertexToCellIndices( aLocalVertexIndex );
    }

    uint Surface_Mesh::get_number_of_cells() const
    {
        return static_cast< uint >( mLocalToGlobalCellIndex.size() );
    }

    uint Surface_Mesh::get_number_of_vertices() const
    {
        return static_cast< uint >( mLocalToGlobalVertexIndex.size() );
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates_of_cell( moris_index aLocalCellIndex ) const
    {
        Matrix< DDRMat >           tVertexCoordinates = this->get_vertex_coordinates();
        moris::Cell< moris_index > tVertexIndices     = this->get_vertices_of_cell( aLocalCellIndex );
        size_t const               tDim               = tVertexCoordinates.n_rows();
        size_t const               tNumVertices       = tVertexIndices.size();
        Matrix< DDRMat >           tCellVertexCoordinates{ tDim, tNumVertices };
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCellVertexCoordinates.set_column( i, tVertexCoordinates.get_column( tVertexIndices( i ) ) );
        }
        return tCellVertexCoordinates;
    }

    Matrix< DDRMat > Surface_Mesh::get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const
    {
        Matrix< DDRMat >           tVertexNormals = this->get_vertex_normals();
        moris::Cell< moris_index > tVertexIndices = this->get_vertices_of_cell( aLocalCellIndex );
        size_t const               tDim           = tVertexNormals.n_rows();
        size_t const               tNumVertices   = tVertexIndices.size();
        Matrix< DDRMat >           tCellVertexNormals{ tDim, tNumVertices };
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCellVertexNormals.set_column( i, tVertexNormals.get_column( tVertexIndices( i ) ) );
        }
        return tCellVertexNormals;
    }

    Json Surface_Mesh::to_json() const
    {
        Json tMesh;

        Json tVertexMap;
        for ( moris::size_t i = 0; i < mLocalToGlobalVertexIndex.size(); i++ )
        {
            moris_index const tGlobalIndex = mLocalToGlobalVertexIndex( i );
            moris_id const    tGlobalID    = mIGMesh->get_mtk_vertex( tGlobalIndex ).get_id();
            tVertexMap.add( std::to_string( tGlobalID ), i );
        }
        tMesh.put_child( "vertex_map", tVertexMap );

        Json tCellMap;
        for ( moris::size_t i = 0; i < mLocalToGlobalCellIndex.size(); i++ )
        {
            moris_index const tGlobalIndex = mLocalToGlobalCellIndex( i );
            moris_id const    tGlobalID    = mIGMesh->get_mtk_cell( tGlobalIndex ).get_id();
            tCellMap.add( std::to_string( tGlobalID ), i );
        }
        tMesh.put_child( "cell_map", tCellMap );

        Json tSideSetMap;
        for ( auto const &tSideSet : mSideSets )
        {
            Json tObj;
            tObj.put( "", tSideSet->get_set_name() );
            tSideSetMap.push_back( { "", tObj } );
        }
        tMesh.put_child( "side_sets", tSideSetMap );

        return tMesh;
    }
}    // namespace moris::mtk