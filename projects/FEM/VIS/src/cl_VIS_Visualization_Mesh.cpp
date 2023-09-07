/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_VIS_Visualization_Mesh.cpp
 *
 */

#include "cl_VIS_Visualization_Mesh.hpp"

namespace moris::vis
{
    // ----------------------------------------------------------------------------

    Visualization_Mesh::Visualization_Mesh( bool aOnlyPrimary )
        : mOnlyPrimary( aOnlyPrimary )
    {
    }

    // ----------------------------------------------------------------------------

    Visualization_Mesh::~Visualization_Mesh()
    {
        // delete sets
        for ( auto iSet : mListOfAllSets )
        {
            delete iSet;
        }

        // delete cells
        for ( auto iCell : mCells )
        {
            delete iCell;
        }

        // delete vertices
        for ( auto iVertex : mVertices )
        {
            delete iVertex;
        }

        // delete clusters
        for ( const auto& iClustersOnSet : mClustersOnBlockSets )
        {
            for ( auto iCluster : iClustersOnSet )
            {
                delete iCluster;
            }
        }

        for ( const auto& iClustersOnSet : mClustersOnSideSets )
        {
            for ( auto iCluster : iClustersOnSet )
            {
                delete iCluster;
            }
        }

        for ( const auto& iClustersOnSet : mLeaderSideClusters )
        {
            for ( auto iCluster : iClustersOnSet )
            {
                delete iCluster;
            }
        }

        for ( const auto& iClustersOnSet : mFollowerSideClusters )
        {
            for ( auto iCluster : iClustersOnSet )
            {
                delete iCluster;
            }
        }

        for ( const auto& iClustersOnSet : mClustersOnDoubleSideSets )
        {
            for ( auto iCluster : iClustersOnSet )
            {
                delete iCluster;
            }
        }
    }

    // ----------------------------------------------------------------------------

    moris::Cell< std::string >
    Visualization_Mesh::get_set_names( mtk::EntityRank aSetEntityRank ) const
    {
        if ( aSetEntityRank == mtk::EntityRank::ELEMENT )
        {
            return mBlockSetNames;
        }
        else if ( aSetEntityRank == mtk::EntityRank::EDGE || aSetEntityRank == mtk::EntityRank::FACE )
        {
            // concatenate list of side set and dbl side sets
            Cell< std::string > tAllOutputSideSetNames = mSideSetNames;
            tAllOutputSideSetNames.append( mDoubleSideSetNames );

            // return the list of all side set names to be outputted
            return tAllOutputSideSetNames;
        }
        else if ( aSetEntityRank == mtk::EntityRank::NODE )
        {
            // don't output node sets
            return {};
        }
        else
        {
            MORIS_ERROR( false,
                    "VIS::Visualization_Mesh::get_set_names() - "
                    "Requested entity rank for set not known. Only block sets and side sets exist." );
            return {};
        }
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_num_nodes() const
    {
        return mVertices.size();
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_num_elems() const
    {
        return mCells.size();
    }

    // ----------------------------------------------------------------------------

    moris::Cell< mtk::Cell const * >
    Visualization_Mesh::get_set_cells( std::string aSetName ) const
    {
        // get the set's position in the list of all sets
        uint tSetIndex = mSetNameToIndexMap.find( aSetName );

        // get access to this set
        moris::mtk::Set* tSet = this->get_set_by_index( tSetIndex );
        MORIS_ERROR(
                tSet->get_set_type() == mtk::SetType::BULK,
                "VIS::Visualization_Mesh::get_set_cells() - "
                "'%s' is not a bulk set. Cannot retrieve cells for this set.",
                aSetName.c_str() );

        // initialize list of cells on set to return
        uint                            tNumCellsOnSet = tSet->get_num_cells_on_set( mOnlyPrimary );
        moris::Cell< const mtk::Cell* > tCellsOnSet( tNumCellsOnSet );

        // collect cells on set
        uint tCellIndexOnSet = 0;
        for ( auto iCluster : tSet->get_clusters_on_set() )
        {
            // collect primary cells in cluster
            for ( auto iPrimaryCell : iCluster->get_primary_cells_in_cluster() )
            {
                tCellsOnSet( tCellIndexOnSet++ ) = iPrimaryCell;
            }

            // if overlapping mesh, also collect the void cells in the cluster
            if ( !mOnlyPrimary )
            {
                for ( auto iVoidCell : iCluster->get_void_cells_in_cluster() )
                {
                    tCellsOnSet( tCellIndexOnSet++ ) = iVoidCell;
                }
            }
        }

        MORIS_ASSERT(
                tCellIndexOnSet == tNumCellsOnSet,
                "VIS::Visualization_Mesh::get_set_cells() - "
                "Number of collected and reported cells on Set '%s' don't match. Something went wrong.",
                aSetName.c_str() );

        // return list of all cells on cluster
        return tCellsOnSet;
    }

    // ----------------------------------------------------------------------------

    moris::Cell< moris::mtk::Vertex const * >
    Visualization_Mesh::get_all_vertices() const
    {
        moris::Cell< moris::mtk::Vertex const * > tVerticesCopy( mVertices.size() );
        for ( uint iVertex = 0; iVertex < mVertices.size(); iVertex++ )
        {
            tVerticesCopy( iVertex ) = mVertices( iVertex );
        }

        return tVerticesCopy;
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Visualization_Mesh::get_num_blocks() const
    {
        return mListOfBlocks.size();
    }

    // ----------------------------------------------------------------------------

    moris::uint
    Visualization_Mesh::get_num_sets() const
    {
        return mListOfAllSets.size();
    }

    // ----------------------------------------------------------------------------

    moris::mtk::Set*
    Visualization_Mesh::get_set_by_index( moris::uint aSetIndex ) const
    {
        MORIS_ASSERT( (uint)aSetIndex < mListOfAllSets.size(), "VIS::Visualization_Mesh::get_set_by_index() - Set index out of bounds" );
        return mListOfAllSets( aSetIndex );
    }

    // ----------------------------------------------------------------------------

    moris::mtk::Set*
    Visualization_Mesh::get_set_by_name( std::string aSetLabel ) const
    {
        MORIS_ASSERT(
                mSetNameToIndexMap.key_exists( aSetLabel ),
                "VIS::Visualization_Mesh::get_set_by_name() - Set with name '%s' is not part of VIS mesh.",
                aSetLabel.c_str() );

        moris_index tSetIndex = mSetNameToIndexMap.find( aSetLabel );
        return mListOfAllSets( tSetIndex );
    }

    // ----------------------------------------------------------------------------

    moris_index
    Visualization_Mesh::get_set_index_by_name( const std::string& aSetLabel ) const
    {
        return mSetNameToIndexMap.find( aSetLabel );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Visualization_Mesh::get_element_indices_in_block_set( uint aSetIndex )
    {
        return mListOfBlocks( aSetIndex )->get_cell_inds_on_block( false );
    }

    // ----------------------------------------------------------------------------

    mtk::MeshType
    Visualization_Mesh::get_mesh_type() const
    {
        return mtk::MeshType::VIS;
    }

    // ----------------------------------------------------------------------------

    Matrix< IdMat >
    Visualization_Mesh::get_communication_table() const
    {
        MORIS_ERROR( false, "VIS::Visualization_Mesh::get_communication_table(), not implemented for visualization mesh" );
        return Matrix< IdMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_spatial_dim() const
    {
        return this->get_set_by_index( 0 )->get_spatial_dim();
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_num_entities(
            mtk::EntityRank aEntityRank,
            moris_index     aDiscretizationIndex ) const
    {
        MORIS_ERROR( false, "VIS::Visualization_Mesh::get_num_entities() - not implemented for visualization mesh" );
        return 0;
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Visualization_Mesh::get_nodes_connected_to_element_loc_inds( moris_index aElementIndex ) const
    {
        return mCells( aElementIndex )->get_vertex_inds();
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Visualization_Mesh::get_entity_connected_to_entity_loc_inds(
            moris_index     aEntityIndex,
            mtk::EntityRank aInputEntityRank,
            mtk::EntityRank aOutputEntityRank,
            moris_index     aDiscretizationIndex ) const
    {
        MORIS_ERROR( false, "VIS::Visualization_Mesh::get_entity_connected_to_entity_loc_inds() - not implemented for visualization mesh" );
        return Matrix< IndexMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Visualization_Mesh::get_elements_connected_to_element_and_face_ord_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( false, "VIS::Visualization_Mesh::get_entity_connected_to_entity_loc_inds() - not implemented for visualization mesh" );
        return Matrix< IndexMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------

    Matrix< IndexMat >
    Visualization_Mesh::get_elements_connected_to_element_and_face_ind_loc_inds( moris_index aElementIndex ) const
    {
        MORIS_ERROR( false, "VIS::Visualization_Mesh::get_entity_connected_to_entity_loc_inds() - not implemented for visualization mesh" );
        return Matrix< IndexMat >( 0, 0 );
    }

    // ----------------------------------------------------------------------------

    Matrix< DDRMat >
    Visualization_Mesh::get_node_coordinate( moris_index aNodeIndex ) const
    {
        return mVertices( aNodeIndex )->get_coords();
    }

    // ----------------------------------------------------------------------------

    moris_id
    Visualization_Mesh::get_glb_entity_id_from_entity_loc_index(
            moris_index     aEntityIndex,
            mtk::EntityRank aEntityRank,
            moris_index     aBSplineMeshIndex ) const
    {
        switch ( aEntityRank )
        {
            case mtk::EntityRank::NODE:
                return mVertices( aEntityIndex )->get_id();
            case mtk::EntityRank::ELEMENT:
                return mCells( aEntityIndex )->get_id();
            default:
                MORIS_ERROR( false, "VIS::Visualization_Mesh::get_glb_entity_id_from_entity_loc_index() - Unknown entity rank." );
                return 0;
        }
    }

    // ----------------------------------------------------------------------------

    moris_index
    Visualization_Mesh::get_loc_entity_ind_from_entity_glb_id(
            moris_id        aEntityId,
            mtk::EntityRank aEntityRank,
            moris_index     aDiscretizationIndex ) const
    {
        MORIS_ERROR( false, "VIS mesh does not yet implement get_loc_entity_ind_from_entity_glb_id()." );
        return 0;
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_node_owner( moris_index aNodeIndex ) const
    {
        return par_rank();
    }

    // ----------------------------------------------------------------------------

    uint
    Visualization_Mesh::get_element_owner( moris_index aElementIndex ) const
    {
        return par_rank();
    }

    // ----------------------------------------------------------------------------

    mtk::CellTopology
    Visualization_Mesh::get_blockset_topology( const std::string& aSetName )
    {
        // TODO: add checks that the set name is actually a blockset + map
        return mListOfBlocks( mSetNameToIndexMap[ aSetName ] )->get_cell_topology();
    }

    // ----------------------------------------------------------------------------

    mtk::CellShape
    Visualization_Mesh::get_IG_blockset_shape( const std::string& aSetName )
    {
        // TODO: add checks that the set name is actually a blockset + map
        return mListOfBlocks( mSetNameToIndexMap[ aSetName ] )->get_IG_cell_shape();
    }

    // ----------------------------------------------------------------------------

    mtk::CellShape
    Visualization_Mesh::get_IP_blockset_shape( const std::string& aSetName )
    {
        // TODO: add checks that the set name is actually a blockset + map
        return mListOfBlocks( mSetNameToIndexMap[ aSetName ] )->get_IP_cell_shape();
    }

    // ----------------------------------------------------------------------------

    void
    Visualization_Mesh::get_sideset_elems_loc_inds_and_ords(
            const std::string&  aSetName,
            Matrix< IndexMat >& aElemIndices,
            Matrix< IndexMat >& aSidesetOrdinals ) const
    {
        // copy string such that we can modify it
        std::string tSetNameCopy = aSetName;

        // get the set
        moris_index tSetIndex = mSetNameToIndexMap.find( tSetNameCopy );
        MORIS_ASSERT(
                (uint)tSetIndex < mListOfAllSets.size(),
                "VIS::Visualization_Mesh::get_sideset_elems_loc_inds_and_ords() - Set index out of bounds" );
        const mtk::Set* tSet = mListOfAllSets( tSetIndex );

        // access the list of side clusters on the current set
        Cell< mtk::Cluster const * > tClustersOnSideSet( 0 );

        // depending on the set type, get the list of clusters for the given set
        switch ( tSet->get_set_type() )
        {
            // bulk sets are not valid for this function
            case mtk::SetType::BULK:
            {
                MORIS_ERROR( false,
                        "VIS::Visualization_Mesh::get_sideset_elems_loc_inds_and_ords() - "
                        "Set '%s' is a block set. This function cannot be applied to block sets.",
                        aSetName.c_str() );
                break;
            }

            // side set
            case mtk::SetType::SIDESET:
            {
                // get the index of the set in the list of side sets
                moris_index tSideSetIndex = mSetLocalToTypeIndex( tSetIndex );

                // get access to the correct list of clusters for output
                tClustersOnSideSet = mClustersOnSideSets( tSideSetIndex );

                // stop switch
                break;
            }

            // double side set
            case mtk::SetType::DOUBLE_SIDED_SIDESET:
            {
                // get the index of the set in the list of double side sets
                moris_index tDblSideSetIndex = mSetLocalToTypeIndex( tSetIndex );

                // use the leader side sets for printing
                tClustersOnSideSet = mLeaderSideClusters( tDblSideSetIndex );

                // stop switch
                break;
            }

            // unknown set type - error
            default:
            {
                MORIS_ERROR( false,
                        "VIS::Visualization_Mesh::get_sideset_elems_loc_inds_and_ords() - "
                        "Set '%s' has unknown set type.",
                        aSetName.c_str() );
            }

        }    // end switch: get clusters depending on set type

        // count number of facets making up the side set
        uint tNumFacets = 0;
        for ( auto iCluster : tClustersOnSideSet )
        {
            tNumFacets += iCluster->get_num_primary_cells();
        }

        // initialize output information
        aElemIndices.set_size( tNumFacets, 1, -1 );
        aSidesetOrdinals.set_size( tNumFacets, 1, -1 );

        // go over clusters and collect all IG cells and side ordinals that make up the side set
        uint tFacetIndex = 0;
        for ( auto iCluster : tClustersOnSideSet )
        {
            // collect the IG elements and ordinals on the current cluster
            Matrix< IndexMat > tIgCellIndicesInCluster = iCluster->get_primary_cell_indices_in_cluster();
            Matrix< IndexMat > tSideOrdinalsInCluster  = iCluster->get_cell_side_ordinals();

            // copy cells and ordinals of the facets into the arrays
            for ( uint iFacetOnCluster = 0; iFacetOnCluster < tSideOrdinalsInCluster.numel(); iFacetOnCluster++ )
            {
                aElemIndices( tFacetIndex )     = tIgCellIndicesInCluster( iFacetOnCluster );
                aSidesetOrdinals( tFacetIndex ) = tSideOrdinalsInCluster( iFacetOnCluster );
                tFacetIndex++;
            }
        }

    }

    void
    Visualization_Mesh::collect_all_sets()
    {
        // collect all mesh sets in one sequential list
        mListOfAllSets.append( mListOfBlocks );
        mListOfAllSets.append( mListOfSideSets );
        mListOfAllSets.append( mListOfDoubleSideSets );

        // construct map linking set name to position in list of all sets
        for ( uint iSet = 0; iSet < mListOfAllSets.size(); iSet++ )
        {
            // get the current set's name
            const std::string tSetName = mListOfAllSets( iSet )->get_set_name();

            // populate the map
            mListOfAllSets( iSet )->set_set_index( iSet );
            mSetNameToIndexMap[ tSetName ] = iSet;
        }
    }

    // ----------------------------------------------------------------------------

    void
    Visualization_Mesh::finalize()
    {
        // collect all types of sets in one sequential list and construct map
        this->collect_all_sets();

        // construct a map that links set names to types and position within that type's list
        mSetLocalToTypeIndex.resize( mListOfAllSets.size() );
        for ( uint iBlockSet = 0; iBlockSet < mBlockSetNames.size(); iBlockSet++ )
        {
            std::string tBlockSetName               = mBlockSetNames( iBlockSet );
            moris_index tGlobalSetIndex             = mSetNameToIndexMap.find( tBlockSetName );
            mSetLocalToTypeIndex( tGlobalSetIndex ) = iBlockSet;
        }
        for ( uint iSideSet = 0; iSideSet < mSideSetNames.size(); iSideSet++ )
        {
            std::string tSideSetName                = mSideSetNames( iSideSet );
            moris_index tGlobalSetIndex             = mSetNameToIndexMap.find( tSideSetName );
            mSetLocalToTypeIndex( tGlobalSetIndex ) = iSideSet;
        }
        for ( uint iDblSideSet = 0; iDblSideSet < mDoubleSideSetNames.size(); iDblSideSet++ )
        {
            std::string tDblSideSetName             = mDoubleSideSetNames( iDblSideSet );
            moris_index tGlobalSetIndex             = mSetNameToIndexMap.find( tDblSideSetName );
            mSetLocalToTypeIndex( tGlobalSetIndex ) = iDblSideSet;
        }

        // mark this mesh as finalized
        mMeshIsFinalized = true;
    }
}
