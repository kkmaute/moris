/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersec_Mesh_Data.cpp
 *
 */

#include "cl_MTK_Intersec_Mesh_Data.hpp"

namespace moris
{
    namespace mtk
    {
        //---------------------------------------------------------------------

        Mesh_Intersection_Data::Mesh_Intersection_Data()
                : mFirstExtEntityInds( (moris::moris_index)EntityRank::UNDEFINED, std::numeric_limits< moris::moris_index >::max() )
                , mLocalToGlobalExtNodes( 1, 1 )
                , mFirstAvailableInds( (moris::moris_index)EntityRank::UNDEFINED, std::numeric_limits< moris::moris_index >::max() )
                , mExternalEntities( (moris::moris_index)EntityRank::UNDEFINED, Vector< mesh::Entity >( 0 ) )
        {
        }

        //---------------------------------------------------------------------

        void
        Mesh_Intersection_Data::set_up_external_entity_data( moris::mtk::Mesh* aMeshData )
        {
            // Initialize the member variables
            mFirstAvailableIds.resize( 7, std::numeric_limits< moris::moris_index >::max() );
            mFirstExtEntityInds.resize( 7, std::numeric_limits< moris::moris_index >::max() );
            mFirstAvailableInds.resize( 7, std::numeric_limits< moris::moris_index >::max() );

            int tProcessorRank = moris::par_rank();

            moris::moris_id tFirstNode = aMeshData->get_max_entity_id( EntityRank::NODE ) + 1;
            // moris::moris_id tFirstEdge = aMeshData->get_max_entity_id(EntityRank::EDGE);
            // moris::moris_id tFirstFace = aMeshData->get_max_entity_id(EntityRank::FACE);
            moris::moris_id tFirstElem = aMeshData->get_max_entity_id( EntityRank::ELEMENT ) + 1;

            if ( tProcessorRank == 0 )
            {
                // Processor 1 (rank 0) is responsible for first available Ids
                mFirstAvailableIds( 0 ) = tFirstNode;
                mFirstAvailableIds( 1 ) = MORIS_ID_MAX;
                mFirstAvailableIds( 2 ) = MORIS_ID_MAX;
                mFirstAvailableIds( 3 ) = tFirstElem;
            }

            mFirstExtEntityInds( 0 ) = aMeshData->get_num_entities( mtk::EntityRank::NODE );
            // mFirstExtEntityInds(1) = aMeshData->get_num_entities(mtk::EntityRank::EDGE);
            //  mFirstExtEntityInds(2) = aMeshData->get_num_entities(mtk::EntityRank::FACE);
            mFirstExtEntityInds( 3 ) = aMeshData->get_num_entities( mtk::EntityRank::ELEMENT );

            mFirstAvailableInds( 0 ) = mFirstExtEntityInds( 0 );
            // mFirstAvailableInds(1) = mFirstExtEntityInds(1);
            // mFirstAvailableInds(2) = mFirstExtEntityInds(2);
            mFirstAvailableInds( 3 ) = mFirstExtEntityInds( 3 );
        }

        //---------------------------------------------------------------------

        void
        Mesh_Intersection_Data::update_first_available_index_external_data( moris::moris_index aNewFirstAvailableIndex, enum EntityRank aEntityRank )
        {
            mFirstAvailableInds( (moris::moris_index)aEntityRank ) = aNewFirstAvailableIndex;

            // needs to be deleted later
            moris::moris_index tEntRankInd  = (moris::moris_index)aEntityRank;
            moris::size_t      tInitialSize = mExternalEntities( tEntRankInd ).size();
            moris::size_t      tAddSize     = 1;
            mExternalEntities( tEntRankInd ).resize( tInitialSize + tAddSize, mesh::Entity() );
        }

        //---------------------------------------------------------------------

        void
        Mesh_Intersection_Data::batch_create_new_nodes_external_data(
                Vector< moris_index > const &                    aNewNodeIds,
                Vector< moris_index > const &                    aNewNodeIndices,
                Vector< moris_index > const &                    aNewNodeOwners,
                Vector< moris::Matrix< moris::DDRMat > > const & aNewNodeCoordinates )
        {
            moris::moris_index tEntRankInd  = (moris::moris_index)EntityRank::NODE;
            moris::size_t      tAddSize     = aNewNodeIds.size();
            moris::size_t      tInitialSize = mExternalEntities( tEntRankInd ).size();

            // Initialize
            moris::size_t      j      = 0;
            moris::moris_index tInd   = MORIS_INDEX_MAX;
            moris::moris_id    tId    = MORIS_ID_MAX;
            moris::moris_id    tOwner = MORIS_ID_MAX;

            // Resize
            mExternalEntities( tEntRankInd ).resize( tInitialSize + tAddSize, mesh::Entity() );
            mLocalToGlobalExtNodes.resize( 1, tInitialSize + tAddSize );

            for ( moris::size_t i = tInitialSize; i < tAddSize + tInitialSize; i++ )
            {
                // Add information to entities
                tInd   = aNewNodeIndices( j );
                tId    = aNewNodeIds( j );
                tOwner = aNewNodeOwners( j );

                mLocalToGlobalExtNodes( tInd - mFirstExtEntityInds( 0 ) ) = tId;
                mExternalEntities( tEntRankInd )( tInd - mFirstExtEntityInds( 0 ) ).set_entity_identifiers( tId, tInd, tOwner, mtk::EntityRank::NODE );
                mExternalEntities( tEntRankInd )( tInd - mFirstExtEntityInds( 0 ) ).set_entity_coords( aNewNodeCoordinates( j ) );
                j++;
            }
        }

        //---------------------------------------------------------------------

        moris::moris_id
        Mesh_Intersection_Data::get_glb_entity_id_from_entity_loc_index_external_data(
                moris::moris_id aEntityIndex,
                enum EntityRank aEntityRank ) const
        {
            moris::size_t tEntityRankIndex = (moris::size_t)aEntityRank;

            moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds( tEntityRankIndex );

            return mExternalEntities( tEntityRankIndex )( tExternalIndex ).get_entity_glb_id();
        }

        //---------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > const &
        Mesh_Intersection_Data::get_selected_node_coordinates_loc_inds_external_data( moris::moris_index aEntityIndex ) const
        {
            moris::size_t      tEntityRankIndex = (moris::size_t)EntityRank::NODE;
            moris::moris_index tExternalIndex   = aEntityIndex - mFirstExtEntityInds( tEntityRankIndex );

            return mExternalEntities( tEntityRankIndex )( tExternalIndex ).get_entity_coords();
        }

        //---------------------------------------------------------------------

        void
        Mesh_Intersection_Data::get_all_node_coordinates_loc_inds_external_data( moris::moris_index aStartingIndex,
                moris::Matrix< moris::DDRMat >&                                                     aCoordinates ) const
        {
            moris::size_t tEntityRankIndex = (moris::size_t)EntityRank::NODE;
            moris::size_t tNumNodes        = this->get_num_entities_external_data( EntityRank::NODE );

            for ( moris::moris_index i = 0; i < (moris::moris_index)tNumNodes; i++ )
            {
                const moris::Matrix< moris::DDRMat >& tCoordinateRow = mExternalEntities( tEntityRankIndex )( i ).get_entity_coords();
                aCoordinates.set_row( aStartingIndex, tCoordinateRow );
                aStartingIndex++;
            }
        }

        //---------------------------------------------------------------------

        moris::moris_id
        Mesh_Intersection_Data::allocate_entity_ids_external_entity_data(
                moris::size_t   aNumIdstoAllocate,
                enum EntityRank aEntityRank ) const
        {
            int tProcRank = par_rank();
            int tProcSize = par_size();

            // size_t is defined as uint here because of aNumRequested
            // Initialize gathered information outputs (information which will be scattered across processors)
            Vector< moris::moris_id > aGatheredInfo;
            Vector< moris::moris_id > tFirstId( 1 );
            Vector< moris::moris_id > tNumIdsRequested( 1 );

            tNumIdsRequested( 0 ) = (moris::moris_id)aNumIdstoAllocate;

            moris::gather_vector( tNumIdsRequested, aGatheredInfo );

            Vector< moris::moris_id > tProcFirstID( tProcSize );

            moris::moris_id tEntityRankIndex = (moris::moris_id)aEntityRank;

            if ( tProcRank == 0 )
            {
                // Loop over entities print the number of entities requested by each processor
                for ( int iProc = 0; iProc < tProcSize; ++iProc )
                {
                    // Give each processor their desired amount of IDs
                    tProcFirstID( iProc ) = mFirstAvailableIds( tEntityRankIndex );

                    // Increment the first available node ID
                    mFirstAvailableIds( tEntityRankIndex ) = mFirstAvailableIds( tEntityRankIndex ) + aGatheredInfo( iProc );
                }
            }

            moris::scatter_vector( tProcFirstID, tFirstId );

            return tFirstId( 0 );
        }

        //---------------------------------------------------------------------

        //---------------------------------------------------------------------
    }    // namespace mtk
}    // namespace moris
