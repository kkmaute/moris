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

#include "cl_MTK_Cell_Info.hpp"
#include "cl_MTK_Interpolation_Function.hpp"
#include "cl_MTK_Interpolation_Rule.hpp"

namespace moris::mtk
{
    Integration_Surface_Mesh::Integration_Surface_Mesh(
            Integration_Surface_Mesh_Data const &aData )
            : Surface_Mesh(
                      aData.get_vertex_coordinates(),
                      aData.get_cell_to_vertex_indices(),
                      1e-9 )
            , mData( aData )
            , mFieldInterpRule( get_ip_field_interpolation_rule( aData ) )
            , mFieldSpaceInterpolator( mFieldInterpRule )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

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

    Vector< Vector< moris_index > > Integration_Surface_Mesh::get_vertex_neighbors() const
    {
        return mData.mVertexNeighbors;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Integration_Surface_Mesh::get_vertex_neighbors( moris_index aLocalVertexIndex ) const
    {
        MORIS_ASSERT( aLocalVertexIndex < static_cast< moris_index >( mData.mVertexNeighbors.size() ), "Vertex index out of bounds" );
        return mData.mVertexNeighbors( aLocalVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Integration_Surface_Mesh::get_global_vertex_index( moris_index aLocalVertexIndex ) const
    {
        return mData.mLocalToGlobalVertexIndex( aLocalVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Integration_Surface_Mesh::get_global_cell_index( moris_index aLocalCellIndex ) const
    {
        return mData.mLocalToGlobalCellIndex( aLocalCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Integration_Surface_Mesh::get_local_vertex_index( moris_index aGlobalVertexIndex ) const
    {
        return mData.mGlobalToLocalVertexIndex.at( aGlobalVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Integration_Surface_Mesh::get_local_cell_index( moris_index aGlobalCellIndex ) const
    {
        return mData.mGlobalToLocalCellIndex.at( aGlobalCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    const std::map< moris_index, Matrix< DDRMat > > &Integration_Surface_Mesh::get_ip_element_displacement() const
    {
        return mIPElementDisplacements;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Integration_Surface_Mesh::set_ip_element_displacement( moris_index aIPCellIndex, const Matrix< DDRMat > &aIPElementDisplacements )
    {
        // Store the displacement of the IP element corresponding to the given facet index
        mIPElementDisplacements[ aIPCellIndex ] = aIPElementDisplacements;

        // -------------------------------------------------------------------
        // Update facet displacement given the IP element displacement
        // -------------------------------------------------------------------

        // Loop through all facets
        for ( uint iF = 0; iF < this->get_number_of_facets(); iF++ )
        {
            if ( mData.mIPElementFacetIndex( iF ) == aIPCellIndex )
            {
                Vector< moris_index > tVertexIndices = this->get_facets_vertex_indices( iF );
                size_t const          tNumVertices   = tVertexIndices.size();

                // Get information about the IP cell and cluster for this facet
                const mtk::Cluster *tCluster                 = mData.mFacetClusters( iF );
                uint                tLeaderClusterLocalIndex = mData.mIPClusterLocalIndex( iF );

                if ( tCluster )
                {
                    // Now use this index for local coordinates
                    Matrix< DDRMat > tTargetLocalCoordinates = tCluster->get_cell_local_coords_on_side_wrt_interp_cell( tLeaderClusterLocalIndex );

                    // Get the interpolation cell from the cluster
                    const mtk::Cell &tIPElement = tCluster->get_interpolation_cell();

                    // Get the local coordinates of the vertices of the interpolation cell
                    Matrix< DDRMat >      tIPElementVertices;
                    const mtk::Cell_Info *tIPInfo = tIPElement.get_cell_info();
                    tIPInfo->get_loc_coords_of_cell( tIPElementVertices );

                    // Go through map to find the displacement for this IP element
                    moris_index tIPElementIndex = tIPElement.get_index();
                    auto        tIt             = mIPElementDisplacements.find( tIPElementIndex );
                    if ( tIt != mIPElementDisplacements.end() )
                    {
                        // Get the displacement for this IP element
                        Matrix< DDRMat > tDispForInterp = tIt->second;

                        // Set field space interpolator coefficients
                        mFieldSpaceInterpolator.set_space_coeff( tDispForInterp );
                        mFieldSpaceInterpolator.set_space_param_coeff( tIPElementVertices );

                        // Interpolate displacement to each cell vertex
                        for ( size_t iV = 0; iV < tNumVertices; iV++ )
                        {
                            // Get the local coordinates of the vertex in the parametric space of the IP element
                            Matrix< DDRMat > tVertexParamCoord = trans( tTargetLocalCoordinates.get_row( iV ) );
                            mFieldSpaceInterpolator.set_space( tVertexParamCoord );

                            // Interpolate the displacement for this vertex and set it in the surface mesh
                            Matrix< DDRMat > tVertexDisp = mFieldSpaceInterpolator.valx();    // (dim x 1)
                            this->set_vertex_displacement( tVertexIndices( iV ), tVertexDisp );
                        }
                    }
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    moris_index Integration_Surface_Mesh::get_cluster_of_cell( moris_index aLocalCellIndex ) const
    {
        return mData.mCellToClusterIndices( aLocalCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Integration_Surface_Mesh::get_vertices_of_cell( moris_index aLocalCellIndex ) const
    {
        return mData.mFacetToVertexIndices( aLocalCellIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< moris_index > Integration_Surface_Mesh::get_cells_of_vertex( moris_index aLocalVertexIndex ) const
    {
        return mData.mVertexToCellIndices( aLocalVertexIndex );
    }

    //--------------------------------------------------------------------------------------------------------------

    uint Integration_Surface_Mesh::get_spatial_dimension() const
    {
        return mData.mIGMesh->get_spatial_dim();
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Integration_Surface_Mesh::get_vertex_normals_of_cell( moris_index aLocalCellIndex ) const
    {
        Matrix< DDRMat >      tVertexNormals = this->get_vertex_normals();
        Vector< moris_index > tVertexIndices = this->get_vertices_of_cell( aLocalCellIndex );
        size_t const          tDim           = this->get_spatial_dimension();
        size_t const          tNumVertices   = tVertexIndices.size();
        Matrix< DDRMat >      tCellVertexNormals{ tDim, tNumVertices };
        for ( moris::size_t i = 0; i < tNumVertices; i++ )
        {
            tCellVertexNormals.set_column( i, tVertexNormals.get_column( tVertexIndices( i ) ) );
        }
        return tCellVertexNormals;
    }

    //--------------------------------------------------------------------------------------------------------------

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

    //--------------------------------------------------------------------------------------------------------------

    Interpolation_Rule Integration_Surface_Mesh::get_ip_field_interpolation_rule(
            Integration_Surface_Mesh_Data const &aData )
    {
        for ( uint iF = 0; iF < aData.get_number_of_facets(); iF++ )
        {
            const Cluster *tCluster = aData.mFacetClusters( iF );

            if ( tCluster != nullptr )
            {
                const Cell &tIPElement =
                        tCluster->get_interpolation_cell();

                return Interpolation_Rule(
                        tIPElement.get_geometry_type(),
                        Interpolation_Type::LAGRANGE,
                        tIPElement.get_cell_info()->get_cell_interpolation_order(),
                        Interpolation_Type::UNDEFINED,
                        Interpolation_Order::UNDEFINED );
            }
        }

        MORIS_ERROR( false, "No facet cluster found for entire surface mesh, so field interpolation rule cannot be found" );

        // Unreachable
        return Interpolation_Rule(
                Geometry_Type::UNDEFINED,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED,
                Interpolation_Type::UNDEFINED,
                Interpolation_Order::UNDEFINED );
    }
}    // namespace moris::mtk
