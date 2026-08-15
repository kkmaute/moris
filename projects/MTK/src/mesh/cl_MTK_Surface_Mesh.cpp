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
#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Interpolation_Function.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Space_Interpolator.hpp"

namespace moris::mtk
{
    Surface_Mesh::Surface_Mesh(
            Integration_Mesh const      *aIGMesh,
            const Vector< std::string > &aSideSetNames )
            : Surface_Mesh( dynamic_cast< Integration_Mesh_DataBase_IG const * >( aIGMesh ), aSideSetNames )
    {
    }

    Surface_Mesh::Surface_Mesh(
            Integration_Mesh const     *aIGMesh,
            Vector< Side_Set const * > &aSideSets )
            : mIGMesh( dynamic_cast< Integration_Mesh_DataBase_IG const * >( aIGMesh ) )
    {
        this->initialize_from_side_sets( aSideSets );
    }

    Surface_Mesh::Surface_Mesh(
            Integration_Mesh_DataBase_IG const *aIGMesh,
            const Vector< std::string >        &aSideSetNames )
    {
        Vector< Side_Set const * > tSideSets;

        auto tSideSetFromName = [ &aIGMesh ]( std::string const &aSideSetName ) {
            return dynamic_cast< Side_Set * >( aIGMesh->get_set_by_name( aSideSetName ) );
        };

        std::transform( aSideSetNames.begin(), aSideSetNames.end(), std::back_inserter( tSideSets ), tSideSetFromName );

        initialize_from_side_sets( tSideSets );
    }

    void Surface_Mesh::initialize_from_side_sets( Vector< Side_Set const * > const &aSideSets )
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
        this->initialize_vertex_coordinates();
        this->initialize_facet_normals();
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    void Surface_Mesh::initialize_side_set(
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

    void Surface_Mesh::initialize_cluster(
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

    void Surface_Mesh::initialize_cell(
            map< moris_index, Vector< moris_index > > &aTmpNeighborMap,
            const Cell                                *aCell,
            int                                        aCellOrdinal,
            moris_index                                aClusterIndex )
    {
        if ( mGlobalToLocalCellIndex.count( aCell->get_index() ) != 0 )
        {
            // if the cell has already been added to the surface mesh, we do not add it again as it would break the mGlobalToLocalCellIndex map
            // FIXME: not correct way to handle; a cell can have two sides on the surface mesh
            //        if ( mCellToClusterIndices( mGlobalToLocalCellIndex[ aCell->get_index() ] ) != aClusterIndex )
            //        {
            //            MORIS_ERROR( false, "Cell already exists in surface mesh with different cluster index" );
            //        }
            //            fprintf( stderr, "Cell %d already exists in surface mesh, side ordinal %d, skipping initialization",    //
            //                    aCell->get_index(),
            //                    aCellOrdinal );
            return;
        }
        // MORIS_ASSERT( mGlobalToLocalCellIndex.count( aCell->get_index() ) == 0, "Cell added twice to surface mesh" );

        // fprintf( stderr, "Adding cell %d to surface mesh with cluster index %d and side ordinal %d\n", aCell->get_index(), aClusterIndex, aCellOrdinal );

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

    void Surface_Mesh::initialize_vertex(
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

    void Surface_Mesh::initialize_neighbors( map< moris_index, Vector< moris_index > > &aTmpNeighborMap )
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
        uint const tDim         = this->get_spatial_dimension();
        mVertexCoordinates.resize( tDim, tNumVertices );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            mVertexCoordinates.set_column( i, mIGMesh->get_node_coordinate( mLocalToGlobalVertexIndex( i ) ) );
        }
    }

    void Surface_Mesh::initialize_facet_normals()
    {
        // Compute facet normals using deformed vertex coordinates.
        auto const tNumCells = static_cast< moris::size_t >( mLocalToGlobalCellIndex.size() );
        uint const tDim      = this->get_spatial_dimension();
        mFacetNormals.resize( tDim, tNumCells );
        for ( moris::size_t i = 0; i < tNumCells; i++ )
        {
            Matrix< DDRMat > tCoords = get_vertex_coordinates_of_cell( i );    // use deformed coordinates
            Matrix< DDRMat > tNormal( tDim, 1 );

            if ( tDim == 2 )
            {
                // 2D: rotate tangent vector
                tNormal( 0 ) = tCoords( 1, 1 ) - tCoords( 1, 0 );
                tNormal( 1 ) = tCoords( 0, 0 ) - tCoords( 0, 1 );
            }
            else    // 3D
            {
                // 3D: cross product of two edge vectors
                Matrix< DDRMat > tEdge1 = tCoords.get_column( 1 ) - tCoords.get_column( 0 );
                Matrix< DDRMat > tEdge2 = tCoords.get_column( 2 ) - tCoords.get_column( 0 );
                tNormal( 0 )            = tEdge1( 1 ) * tEdge2( 2 ) - tEdge1( 2 ) * tEdge2( 1 );
                tNormal( 1 )            = tEdge1( 2 ) * tEdge2( 0 ) - tEdge1( 0 ) * tEdge2( 2 );
                tNormal( 2 )            = tEdge1( 0 ) * tEdge2( 1 ) - tEdge1( 1 ) * tEdge2( 0 );
            }
            tNormal = tNormal / norm( tNormal );
            mFacetNormals.set_column( i, tNormal );
        }
    }

    void Surface_Mesh::initialize_facet_measure()
    {
        // Compute facet measures using deformed vertex coordinates.
        auto const tNumCells = static_cast< moris::size_t >( mLocalToGlobalCellIndex.size() );
        uint const tDim      = this->get_spatial_dimension();
        mFacetMeasure.resize( tNumCells, 1 );
        for ( moris::size_t i = 0; i < tNumCells; i++ )
        {
            Matrix< DDRMat > tCoords = get_vertex_coordinates_of_cell( i );
            if ( tDim == 2 )
            {
                // 2D: length of line segment
                mFacetMeasure( i ) = norm( tCoords.get_column( 1 ) - tCoords.get_column( 0 ) );
            }
            else    // 3D
            {
                // 3D: area using cross product (for triangles/quads this is approximate)
                Matrix< DDRMat > tEdge1 = tCoords.get_column( 1 ) - tCoords.get_column( 0 );
                Matrix< DDRMat > tEdge2 = tCoords.get_column( 2 ) - tCoords.get_column( 0 );
                Matrix< DDRMat > tCross( 3, 1 );
                tCross( 0 )        = tEdge1( 1 ) * tEdge2( 2 ) - tEdge1( 2 ) * tEdge2( 1 );
                tCross( 1 )        = tEdge1( 2 ) * tEdge2( 0 ) - tEdge1( 0 ) * tEdge2( 2 );
                tCross( 2 )        = tEdge1( 0 ) * tEdge2( 1 ) - tEdge1( 1 ) * tEdge2( 0 );
                mFacetMeasure( i ) = 0.5 * norm( tCross );    // Triangle area
                // For quads with 4 vertices, would need to sum two triangles
            }
        }
    }

    void Surface_Mesh::initialize_vertex_normals()
    {
        auto const tNumVertices = static_cast< moris::size_t >( mLocalToGlobalVertexIndex.size() );
        uint const tDim         = this->get_spatial_dimension();
        mVertexNormals.resize( tDim, tNumVertices );
        auto tNormal = Matrix< DDRMat >( tDim, 1 );
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            Vector< moris_index > tVertexCellNeighbors = mVertexToCellIndices( i );
            auto const            tNumNeighbors        = static_cast< moris::size_t >( tVertexCellNeighbors.size() );
            tNormal.fill( 0.0 );
            // compute the normal as the weighted average of the facet normals of the neighboring cells
            for ( moris::size_t j = 0; j < tNumNeighbors; j++ )
            {
                int const tCellIndex = tVertexCellNeighbors( j );
                // Use deformed facet measure and normals
                tNormal += mFacetNormals.get_column( tCellIndex ) * mFacetMeasure( tCellIndex );
            }
            mVertexNormals.set_column( i, tNormal / norm( tNormal ) );
        }
    }

    // void Surface_Mesh::interpolate_facet_vertex_displacements()
    // {
    //     for (const auto& [tIPElementIndex, tDisp] : mFacetDisplacements)
    //     {
    //         // Find the cluster for this interpolation element
    //         const mtk::Cluster* tCluster = nullptr;
    //         for (const Side_Set* tSideSet : mSideSets)
    //         {
    //             for (uint iCluster = 0; iCluster < tSideSet->get_num_clusters_on_set(); ++iCluster)
    //             {
    //                 const mtk::Cluster* candidateCluster = tSideSet->get_clusters_by_index(iCluster);
    //                 const mtk::Cell& tIPElement = candidateCluster->get_interpolation_cell();
    //                 if (tIPElement.get_index() == tIPElementIndex)
    //                 {
    //                     tCluster = candidateCluster;
    //                     break;
    //                 }
    //             }
    //             if (tCluster) break;
    //         }
    //         if (!tCluster) continue;
    //         // Get local coordinates of facet vertices
    //         auto numPrimaryCells = tCluster->get_num_primary_cells();
    //         for (uint cellIdx = 0; cellIdx < numPrimaryCells; ++cellIdx)
    //         {
    //             Matrix<DDRMat> tTargetLocalCoordinates = tCluster->get_cell_local_coords_on_side_wrt_interp_cell(cellIdx);
    //             size_t tNumVertices = tTargetLocalCoordinates.n_cols();
    //             // Set up cubic interpolation rule
    //             Interpolation_Rule tFieldInterpRule(
    //                 Geometry_Type::LINE,
    //                 Interpolation_Type::LAGRANGE,
    //                 Interpolation_Order::CUBIC,
    //                 Interpolation_Type::UNDEFINED,
    //                 Interpolation_Order::UNDEFINED);
    //             Space_Interpolator tFieldSpaceInterpolator(tFieldInterpRule);
    //             // Ensure tDisp is (dim, numNodes)
    //             Matrix<DDRMat> tDispForInterp = tDisp;
    //             if (tDisp.n_rows() != this->get_spatial_dimension())
    //                 tDispForInterp = trans(tDisp);
    //             tFieldSpaceInterpolator.set_space_coeff(tDispForInterp);
    //             // For each facet vertex, interpolate displacement
    //             for (size_t i = 0; i < tNumVertices; ++i)
    //             {
    //                 Matrix<DDRMat> tVertexParamCoord = tTargetLocalCoordinates.get_column(i);
    //                 tFieldSpaceInterpolator.set_space(tVertexParamCoord);
    //                 Matrix<DDRMat> u_vertex = tFieldSpaceInterpolator.valx();
    //                 // You can now use u_vertex as the interpolated displacement at this facet vertex
    //                 // e.g., store or print as needed
    //             }
    //         }
    //     }
    // }

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates() const
    {
        // Return deformed vertex coordinates when IP element displacements are present.
        Matrix< DDRMat > tVertexCoordinates = mVertexCoordinates;
        // If facet displacements exist, return deformed coordinates for each cell
        if ( !mIPElementDisplacements.empty() )
        {
            Matrix< DDRMat > tDeformedCoordinates = tVertexCoordinates;
            for ( uint localCellIndex = 0; localCellIndex < get_number_of_cells(); ++localCellIndex )
            {
                Matrix< DDRMat >      tCellVertexCoordinates = get_vertex_coordinates_of_cell( localCellIndex );
                Vector< moris_index > tVertexIndices         = get_vertices_of_cell( localCellIndex );
                size_t                tNumVertices           = tVertexIndices.size();
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    tDeformedCoordinates.set_column( tVertexIndices( i ), tCellVertexCoordinates.get_column( i ) );
                }
            }
            return tDeformedCoordinates;
        }
        return tVertexCoordinates;
    }

    Vector< Vector< moris_index > > Surface_Mesh::get_vertex_neighbors() const
    {
        return mVertexNeighbors;
    }

    Vector< moris_index > Surface_Mesh::get_vertex_neighbors( moris_index aLocalVertexIndex ) const
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

    Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates_of_cell( moris_index aLocalCellIndex ) const
    {
        // Compute cell vertex coordinates, applying IP element displacements if available.
        Matrix< DDRMat >      tVertexCoordinates = mVertexCoordinates;
        Vector< moris_index > tVertexIndices     = this->get_vertices_of_cell( aLocalCellIndex );
        size_t const          tDim               = tVertexCoordinates.n_rows();
        size_t const          tNumVertices       = tVertexIndices.size();
        Matrix< DDRMat >      tCellVertexCoordinates{ tDim, tNumVertices };

        // Check if facet displacements exist for this cell
        moris_index tGlobalCellIndex = this->get_global_cell_index( aLocalCellIndex );

        // Use the interpolation element index for facet displacements
        const mtk::Cluster *tCluster = nullptr;
        for ( uint iSideSet = 0; iSideSet < mSideSets.size(); ++iSideSet )
        {
            const Side_Set *tSideSet = mSideSets( iSideSet );
            for ( uint iCluster = 0; iCluster < tSideSet->get_num_clusters_on_set(); ++iCluster )
            {
                const mtk::Cluster *tCandidateCluster = tSideSet->get_clusters_by_index( iCluster );
                auto                tPrimaryCells     = tCandidateCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                for ( uint i = 0; i < tPrimaryCells.size(); ++i )
                {
                    if ( tPrimaryCells( i )->get_index() == tGlobalCellIndex )
                    {
                        tCluster = tCandidateCluster;
                        break;
                    }
                }
                if ( tCluster ) break;
            }
            if ( tCluster ) break;
        }

        if ( tCluster )
        {
            // Find the local index of the cell within the cluster
            moris_index aLeaderClusterLocalIndex = -1;
            auto        numPrimaryCells          = tCluster->get_num_primary_cells();
            for ( uint i = 0; i < numPrimaryCells; ++i )
            {
                Vector< moris::mtk::Cell const * > const &tPrimaryCellsInCluster = tCluster->get_primary_cells_in_cluster( mtk::Leader_Follower::LEADER );
                if ( tPrimaryCellsInCluster( i )->get_index() == tGlobalCellIndex )
                {
                    aLeaderClusterLocalIndex = i;
                    break;
                }
            }

            // Now use this index for local coordinates
            Matrix< DDRMat > tTargetLocalCoordinates = tCluster->get_cell_local_coords_on_side_wrt_interp_cell( aLeaderClusterLocalIndex );
            // Get the interpolation cell from the cluster
            const mtk::Cell &tIPElement = tCluster->get_interpolation_cell();

            Matrix< DDRMat >      tIPElementVertices;
            const mtk::Cell_Info *tIPInfo = tIPElement.get_cell_info();
            tIPInfo->get_loc_coords_of_cell( tIPElementVertices );

            // Get displacement for THIS specific IP element
            moris_index tIPElementIndex = tIPElement.get_index();
            auto        it              = mIPElementDisplacements.find( tIPElementIndex );
            if ( it != mIPElementDisplacements.end() )
            {
                Matrix< DDRMat > tDispForInterp = it->second;
                // Set up space interpolator for displacement using all IP element nodes
                Interpolation_Rule tFieldInterpRule(
                        tIPElement.get_geometry_type(),
                        Interpolation_Type::LAGRANGE,
                        tIPElement.get_cell_info()->get_cell_interpolation_order(),
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
                Space_Interpolator tFieldSpaceInterpolator( tFieldInterpRule );
                tFieldSpaceInterpolator.set_space_coeff( tDispForInterp );
                tFieldSpaceInterpolator.set_space_param_coeff( tIPElementVertices );

                // Interpolate displacement to each cell vertex
                for ( size_t i = 0; i < tNumVertices; i++ )
                {
                    Matrix< DDRMat > tVertexParamCoord = trans( tTargetLocalCoordinates.get_row( i ) );
                    tFieldSpaceInterpolator.set_space( tVertexParamCoord );
                    Matrix< DDRMat > tVertexDisp = tFieldSpaceInterpolator.valx();    // (dim x 1)
                    tCellVertexCoordinates.set_column( i, tVertexCoordinates.get_column( tVertexIndices( i ) ) + trans( tVertexDisp ) );
                }
                return tCellVertexCoordinates;
            }
        }

        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCellVertexCoordinates.set_column( i, tVertexCoordinates.get_column( tVertexIndices( i ) ) );
        }
        return tCellVertexCoordinates;
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
        MORIS_ASSERT( aDisplacements.n_rows() == this->get_spatial_dimension(), "Number of vertices in displacement matrix does not match number of vertices in mesh" );
        MORIS_ASSERT( aDisplacements.n_cols() == mLocalToGlobalVertexIndex.size(), "Number of dimensions in displacement matrix does not match number of dimensions in mesh" );
        mDisplacements = aDisplacements;

        // the displacement on each vertex invalidates the facet and vertex normals as well as the facet measure.
        this->initialize_facet_normals();
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    void Surface_Mesh::set_ip_element_displacement( moris_index aCellIndex, Matrix< DDRMat > const &aIPElementDisplacements )
    {
        mIPElementDisplacements[ aCellIndex ] = aIPElementDisplacements;
    }

    void Surface_Mesh::refresh_derived_quantities()
    {
        this->initialize_facet_normals();
        this->initialize_facet_measure();
        this->initialize_vertex_normals();
    }

    moris_index Surface_Mesh::get_cluster_of_cell( moris_index aLocalCellIndex ) const
    {
        return mCellToClusterIndices( aLocalCellIndex );
    }

    Vector< moris_index > Surface_Mesh::get_vertices_of_cell( moris_index aLocalCellIndex ) const
    {
        return mCellToVertexIndices( aLocalCellIndex );
    }

    Vector< moris_index > Surface_Mesh::get_cells_of_vertex( moris_index aLocalVertexIndex ) const
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

    uint Surface_Mesh::get_spatial_dimension() const
    {
        return mIGMesh->get_spatial_dim();
    }

    void Surface_Mesh::write_to_file( const std::string &aFilePath ) const
    {
        // Open file for writing
        std::ofstream tFile;
        tFile.open( aFilePath );
        tFile << std::fixed << std::setprecision( 8 );

        Matrix< DDRMat > tVertexCoordinates = this->get_vertex_coordinates();

        // Write vertices
        for ( uint iVertex = 0; iVertex < this->get_number_of_vertices(); iVertex++ )
        {
            tFile << "v ";
            for ( uint iDimension = 0; iDimension < this->get_spatial_dimension(); iDimension++ )
            {
                tFile << tVertexCoordinates( iDimension, iVertex ) << " ";
            }
            tFile << "\n";
        }

        for ( uint iFacet = 0; iFacet < this->get_number_of_cells(); iFacet++ )
        {
            tFile << "f ";
            const Vector< moris_index > &tVertexIndices = this->get_vertices_of_cell( iFacet );
            for ( uint iVertexIndex = 0; iVertexIndex < tVertexIndices.size(); iVertexIndex++ )
            {
                tFile << tVertexIndices( iVertexIndex ) + 1 << " ";
            }
            tFile << "\n";
        }

        // close file
        tFile.close();
    }

    // Matrix< DDRMat > Surface_Mesh::get_vertex_coordinates_of_cell( moris_index aLocalCellIndex ) const
    // {
    //     Matrix< DDRMat >      tVertexCoordinates = this->get_vertex_coordinates();
    //     Vector< moris_index > tVertexIndices     = this->get_vertices_of_cell( aLocalCellIndex );
    //     size_t const          tDim               = tVertexCoordinates.n_rows();
    //     size_t const          tNumVertices       = tVertexIndices.size();
    //     Matrix< DDRMat >      tCellVertexCoordinates{ tDim, tNumVertices };
    //     for ( moris::size_t i = 0; i < tNumVertices; i++ )
    //     {
    //         tCellVertexCoordinates.set_column( i, tVertexCoordinates.get_column( tVertexIndices( i ) ) );
    //     }
    //     return tCellVertexCoordinates;
    // }

    Matrix< DDRMat > Surface_Mesh::get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const
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

    // const mtk::Cell &Surface_Mesh::get_global_cell( moris_index aLocalCellIndex ) const
    // {
    //     return mIGMesh->get_mtk_cell( this->get_global_cell_index( aLocalCellIndex ) );
    // }

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

        Json             tCoordinatesJson;
        Matrix< DDRMat > tCoords = get_vertex_coordinates();
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
