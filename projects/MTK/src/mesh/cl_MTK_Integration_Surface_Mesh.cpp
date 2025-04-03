/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Surface_Mesh.cpp
 *
 */

#include "cl_MTK_Integration_Surface_Mesh.hpp"
#include "cl_MTK_Mesh_DataBase_IG.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cell_DataBase.hpp"
#include "fn_trans.hpp"
#include "cl_MTK_Vertex_DataBase.hpp"
#include "cl_Json_Object.hpp"

namespace moris::mtk
{
    Integration_Surface_Mesh::Integration_Surface_Mesh(
            Integration_Surface_Mesh_Data const &aData )
            : Surface_Mesh(
                      aData.get_vertex_coordinates(),
                      aData.get_cell_to_vertex_indices(),
                      1e-9 )
            , mData( aData )
    {
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    Matrix< DDRMat > Integration_Surface_Mesh::initialize_vertex_coordinates( Integration_Mesh const *aIGMesh )
    {
        // uint tNumVertices = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        uint tNumVertices = aIGMesh->get_num_nodes();
        uint tDim         = aIGMesh->get_spatial_dim();

        // Resize the matrix to hold the coordinates
        Matrix< DDRMat > tVertexCoordinates;
        tVertexCoordinates.set_size( tDim, tNumVertices );

        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            // Retrieve the coordinates from mIGMesh and store them in tVertexCoordinates
            tVertexCoordinates.set_column( i, aIGMesh->get_node_coordinate( i ) );
        }

        // Return the coordinates matrix
        return tVertexCoordinates;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Vector< moris_index > > Integration_Surface_Mesh::get_cell_to_vertex_indices(
            Integration_Mesh const           *aIGMesh,
            const Vector< Side_Set const * > &aSideSets ) const
    {
        Vector< Vector< moris_index > > tFacetConnectivity;

        map< moris_index, Vector< moris_index > > tTmpNeighborMap;

        for ( auto const &tSideSet : aSideSets )
        {
            moris_index const tNumClustersOnSet = tSideSet->get_clusters_on_set().size();

            for ( moris_index tClusterIndex = 0; tClusterIndex < tNumClustersOnSet; tClusterIndex++ )
            {
                Cluster const         *tCluster  = tSideSet->get_clusters_by_index( tClusterIndex );
                Vector< const Cell * > tCells    = tCluster->get_primary_cells_in_cluster();
                Matrix< IdMat >        tCellOrds = tCluster->get_cell_side_ordinals();

                for ( uint i = 0; i < tCells.size(); i++ )
                {
                    Cell const *tCurrentCell        = tCells( i );
                    int const   tCurrentCellOrdinal = tCellOrds( i );

                    Vector< moris_index > tCellVertexIndices;

                    Vector< Vertex const * > tSideVertices = tCurrentCell->get_geometric_vertices_on_side_ordinal( tCurrentCellOrdinal );

                    for ( uint j = 0; j < tSideVertices.size(); j++ )
                    {
                        Vertex const *tVertex = tSideVertices( j );
                        tCellVertexIndices.push_back( tVertex->get_index() );
                    }

                    tFacetConnectivity.push_back( tCellVertexIndices );
                }
            }
        }


        return tFacetConnectivity;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< Side_Set const * >
    Integration_Surface_Mesh::obtain_sidesets_from_names( Integration_Mesh_DataBase_IG const *aIGMesh, const Vector< std::string > &aSideSetNames )
    {
        Vector< Side_Set const * > aSideSets;

        auto aSideSetFromName = [ &aIGMesh ]( std::string const &aSideSetName ) -> Side_Set const * {
            return const_cast< Side_Set const * >( dynamic_cast< Side_Set * >( aIGMesh->get_set_by_name( aSideSetName ) ) );
        };

        std::transform( aSideSetNames.begin(), aSideSetNames.end(), std::back_inserter( aSideSets ), aSideSetFromName );

        return aSideSets;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh::initialize_facet_measure()
    {
        auto const       tNumCells          = static_cast< moris::size_t >( mData.mLocalToGlobalCellIndex.size() );
        Matrix< DDRMat > tVertexCoordinates = this->get_all_vertex_coordinates();
        mFacetMeasure.resize( tNumCells, 1 );

        for ( moris::size_t i = 0; i < tNumCells; i++ )
        {
            MORIS_ASSERT( get_spatial_dimension() == 2, "Surface Mesh facet measure only implemented for 2D meshes (Lines)" );
            Vector< moris_index > tVertices = mData.mCellToVertexIndices( i );
            Matrix< DDRMat >      tCoords( 2, 2 );
            // int const   tGlobalCellIndex = mLocalToGlobalCellIndex( i );
            // int         tSideOrdinal     = mCellSideOrdinals( i );
            // Cell const &tCell            = mIGMesh->get_mtk_cell( tGlobalCellIndex );
            // real const  tMeasure         = tCell.compute_cell_side_measure( tSideOrdinal );
            mFacetMeasure( i ) = norm( tVertexCoordinates.get_column( tVertices( 1 ) ) - tVertexCoordinates.get_column( tVertices( 0 ) ) );
        }
    }

    void Integration_Surface_Mesh::initialize_vertex_normals()
    {
        auto const             tNumVertices  = static_cast< moris::size_t >( mData.mLocalToGlobalVertexIndex.size() );
        uint const             tDim          = this->get_spatial_dimension();
        const Matrix< DDRMat > tFacetNormals = this->get_all_facet_normals();
        const Matrix< DDRMat > tFacetMeasure = this->get_facet_measure();
        mVertexNormals.resize( tDim, tNumVertices );

        auto tNormal = Matrix< DDRMat >( tDim, 1 );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            Vector< moris_index > tVertexCellNeighbors = mData.mVertexToCellIndices( i );
            auto const            tNumNeighbors        = static_cast< moris::size_t >( tVertexCellNeighbors.size() );
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

    Vector< Vector< moris_index > > Integration_Surface_Mesh::get_vertex_neighbors() const
    {
        return mData.mVertexNeighbors;
    }

    Vector< moris_index > Integration_Surface_Mesh::get_vertex_neighbors( moris_index aLocalVertexIndex ) const
    {
        MORIS_ASSERT( aLocalVertexIndex < static_cast< moris_index >( mData.mVertexNeighbors.size() ), "Vertex index out of bounds" );
        return mData.mVertexNeighbors( aLocalVertexIndex );
    }

    const Matrix< DDRMat > &Integration_Surface_Mesh::get_facet_measure() const
    {
        return mFacetMeasure;
    }

    Matrix< DDRMat > Integration_Surface_Mesh::get_vertex_normals() const
    {
        return mVertexNormals;
    }

    moris_index Integration_Surface_Mesh::get_global_vertex_index( moris_index aLocalVertexIndex ) const
    {
        return mData.mLocalToGlobalVertexIndex( aLocalVertexIndex );
    }

    moris_index Integration_Surface_Mesh::get_global_cell_index( moris_index aLocalCellIndex ) const
    {
        return mData.mLocalToGlobalCellIndex( aLocalCellIndex );
    }

    moris_index Integration_Surface_Mesh::get_local_vertex_index( moris_index aGlobalVertexIndex ) const
    {
        return mData.mGlobalToLocalVertexIndex.at( aGlobalVertexIndex );
    }

    moris_index Integration_Surface_Mesh::get_local_cell_index( moris_index aGlobalCellIndex ) const
    {
        return mData.mGlobalToLocalCellIndex.at( aGlobalCellIndex );
    }

    void Integration_Surface_Mesh::set_all_displacements( Matrix< DDRMat > const &aDisplacements )
    {
        // Set the displacement of the vertices in the base class
        Surface_Mesh::set_all_displacements( aDisplacements );

        // the displacement on each vertex invalidates the facet and vertex normals as well as the facet measure, so update them.
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    moris_index Integration_Surface_Mesh::get_cluster_of_cell( moris_index aLocalCellIndex ) const
    {
        return mData.mCellToClusterIndices( aLocalCellIndex );
    }

    Vector< moris_index > Integration_Surface_Mesh::get_vertices_of_cell( moris_index aLocalCellIndex ) const
    {
        return mData.mCellToVertexIndices( aLocalCellIndex );
    }

    Vector< moris_index > Integration_Surface_Mesh::get_cells_of_vertex( moris_index aLocalVertexIndex ) const
    {
        return mData.mVertexToCellIndices( aLocalVertexIndex );
    }

    uint Integration_Surface_Mesh::get_spatial_dimension() const
    {
        return mData.mIGMesh->get_spatial_dim();
    }

    Matrix< DDRMat > Integration_Surface_Mesh::get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const
    {
        Matrix< DDRMat >      tVertexNormals = this->get_vertex_normals();
        Vector< moris_index > tVertexIndices = this->get_vertices_of_cell( aLocalCellIndex );
        size_t const          tDim           = tVertexNormals.n_rows();
        size_t const          tNumVertices   = tVertexIndices.size();
        Matrix< DDRMat >      tCellVertexNormals{ tDim, tNumVertices };
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCellVertexNormals.set_column( i, tVertexNormals.get_column( tVertexIndices( i ) ) );
        }
        return tCellVertexNormals;
    }

    Json Integration_Surface_Mesh::to_json() const
    {
        Json tMesh;

        Json tVertexMap;
        for ( moris::size_t i = 0; i < mData.mLocalToGlobalVertexIndex.size(); i++ )
        {
            moris_index const tGlobalIndex = mData.mLocalToGlobalVertexIndex( i );
            moris_id const    tGlobalID    = mData.mIGMesh->get_mtk_vertex( tGlobalIndex ).get_id();
            tVertexMap.add( std::to_string( tGlobalID ), i );
        }
        tMesh.put_child( "vertex_map", tVertexMap );

        Json             tCoordinatesJson;
        Matrix< DDRMat > tCoords = get_all_vertex_coordinates();
        for ( size_t iVertex = 0; iVertex < get_number_of_vertices(); iVertex++ )
        {
            Json tArray;
            for ( size_t iDim = 0; iDim < tCoords.n_rows(); ++iDim )
            {
                Json tObj;
                tObj.put( "", tCoords( iDim, iVertex ) );
                tArray.push_back( { "", tObj } );
            }
            tCoordinatesJson.put_child( std::to_string( iVertex ), tArray );
        }
        tMesh.put_child( "coordinates", tCoordinatesJson );

        Json tCellMap;
        for ( moris::size_t i = 0; i < mData.mLocalToGlobalCellIndex.size(); i++ )
        {
            moris_index const tGlobalIndex = mData.mLocalToGlobalCellIndex( i );
            moris_id const    tGlobalID    = mData.mIGMesh->get_mtk_cell( tGlobalIndex ).get_id();
            tCellMap.add( std::to_string( tGlobalID ), i );
        }
        tMesh.put_child( "cell_map", tCellMap );

        Json tSideSetMap;
        for ( auto const &tSideSet : mData.mSideSets )
        {
            Json tObj;
            tObj.put( "", tSideSet->get_set_name() );
            tSideSetMap.push_back( { "", tObj } );
        }
        tMesh.put_child( "side_sets", tSideSetMap );

        return tMesh;
    }
}    // namespace moris::mtk
