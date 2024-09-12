/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Intersec_Mesh_Data.hpp
 *
 */

#ifndef PROJECTS_MTK_SRC_CL_MTK_INTERSEC_MESH_DATA_HPP_
#define PROJECTS_MTK_SRC_CL_MTK_INTERSEC_MESH_DATA_HPP_

#include <mpi.h>

// TPL: STK Includes

#include "cl_Matrix.hpp"
// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
// XTKL: Containers
#include "cl_Vector.hpp"

#include "assert.hpp"

namespace moris::mtk::mesh
{
    class Entity
    {
      private:
        moris::moris_id                mGlbId;
        moris::moris_id                mLocInd;
        moris::moris_id                mOwningProc;
        EntityRank                     mEntityRank;
        moris::Matrix< moris::DDRMat > mEntityCoordinates;    // If its a node

      public:
        // ----------------------------------------------------------------------------------

        Entity()
                : mGlbId( MORIS_ID_MAX )
                , mLocInd( MORIS_ID_MAX )
                , mOwningProc( MORIS_ID_MAX )
        {
        }

        // ----------------------------------------------------------------------------------

        ~Entity()
        {
        }

        // ----------------------------------------------------------------------------------

        void set_entity_identifiers(
                moris::moris_id aGlbId,
                moris::moris_id aLocInd,
                moris::moris_id aOwnerProc,
                EntityRank      aEntityRank )
        {
            mGlbId      = aGlbId;
            mLocInd     = aLocInd;
            mOwningProc = aOwnerProc;
            mEntityRank = aEntityRank;
        }

        // ----------------------------------------------------------------------------------

        void set_entity_coords( moris::Matrix< moris::DDRMat > const & aCoordinates )
        {
            if ( mEntityRank == mtk::EntityRank::NODE )
            {
                mEntityCoordinates = aCoordinates.copy();
            }
            else
            {
                MORIS_ERROR( false,
                        "Entity::set_entity_coords - Only nodes should have coordinates in this context to avoid duplicate coordinate storage" );
            }
        }

        // ----------------------------------------------------------------------------------

        moris::moris_index get_entity_loc_index() const
        {
            MORIS_ASSERT( mLocInd != std::numeric_limits< moris::moris_id >::max(), "Index has not been set" );

            return mLocInd;
        }

        // ----------------------------------------------------------------------------------

        moris::moris_index get_entity_glb_id() const
        {
            MORIS_ASSERT( mGlbId != std::numeric_limits< moris::moris_id >::max(), "Id has not been set" );
            return mGlbId;
        }

        // ----------------------------------------------------------------------------------

        moris::moris_index
        get_entity_owner() const
        {
            MORIS_ASSERT( mOwningProc != std::numeric_limits< moris::moris_id >::max(), "Owner has not been set" );
            return mOwningProc;
        }

        // ----------------------------------------------------------------------------------

        moris::Matrix< moris::DDRMat > const &
        get_entity_coords() const
        {
            return mEntityCoordinates;
        }

        // ----------------------------------------------------------------------------------

        //                    size_t
        //                    capacity()
        //                    {
        //                        size_t tTotal = 0;
        //                        tTotal += sizeof(mGlbId);
        //                        tTotal += sizeof(mLocInd);
        //                        tTotal += sizeof(mOwningProc);
        //                        tTotal += sizeof(mEntityRank);
        //                        tTotal += mEntityCoordinates.capacity();
        //                        return tTotal;
        //                    }
    };
}    // namespace moris::mtk::mesh

namespace moris::mtk
{
    /**
     * @brief This class is exclusively for the data that is being added to the mesh due to intersection process
     * Namely index, id of the all the vertices and cells created are stored here
     * Also this class is used to allocate new index and ids for the created entities
     */
    class Mesh_Intersection_Data
    {
        // ----------------------------------------------------------------------------------

      private:
        Vector< moris::moris_index > mFirstExtEntityInds;

        // Owned by proc rank 0, other procs UINT_MAX
        // Mutable to preserve const in the allocate entity ids function
        mutable Vector< moris::moris_id > mFirstAvailableIds;

        // Local to Global Node Map
        moris::Matrix< moris::IdMat > mLocalToGlobalExtNodes;

        // Each processor tracks this value
        Vector< moris::moris_index > mFirstAvailableInds;

        // Entity Rank outside, then entity objects inside
        Vector< Vector< mesh::Entity > > mExternalEntities;

      public:
        // ----------------------------------------------------------------------------------
        /*
         *Default Constructor
         */
        Mesh_Intersection_Data();

        // ----------------------------------------------------------------------------------
        /**
         * This function initializes the first id and index of the external data
         * that will be added due to the intersection process
         * @param[ in ] aMeshData An mtk mesh to find id and index
         */
        void
        set_up_external_entity_data( moris::mtk::Mesh* aMeshData );

        // ----------------------------------------------------------------------------------
        /**
         * Get first available index for an entity (node, element, ...)
         * @param[ in ] aEntityRank Entity rank
         */
        moris::moris_index
        get_first_available_index_external_data( enum EntityRank aEntityRank ) const
        {
            return mFirstAvailableInds( (moris::moris_index)aEntityRank );
        }

        // ----------------------------------------------------------------------------------
        /**
         * After an index is allocated to an entity this function is called to update the index count
         * @param[ in ] aNewFirstAvailableIndex The new first available index
         * @param[ in ] aEntityRank Entity rank
         */
        void
        update_first_available_index_external_data( moris::moris_index aNewFirstAvailableIndex, enum EntityRank aEntityRank );

        // ----------------------------------------------------------------------------------
        /*
         * Created a set of new nodes providing the index,id owner and coordinates
         * @param[ in ] aNewNodeIds Node ids
         * @param[ in ] aNewNodeIndices Node indices
         * @param[ in ] aNewNodeOwners Node's owning processor
         * @param[ in ] aNewNodeCoordinates Node coords
         */
        void
        batch_create_new_nodes_external_data(
                Vector< moris_index > const &                    aNewNodeIds,
                Vector< moris_index > const &                    aNewNodeIndices,
                Vector< moris_index > const &                    aNewNodeOwners,
                Vector< moris::Matrix< moris::DDRMat > > const & aNewNodeCoordinates );

        // ----------------------------------------------------------------------------------
        /*
         * Gets number of external entities added
         * @param[ in ] aEntityRank Entity rank
         */
        moris::size_t
        get_num_entities_external_data( enum EntityRank aEntityRank ) const
        {
            return mExternalEntities( (moris::size_t)aEntityRank ).size();
        }

        // ----------------------------------------------------------------------------------
        /**
         * Based on the index of the entity it determines if it is added during intersection process
         * @param[ in ] aEntityIndex index of the entity
         * @param[ in ] aEntityRank rank of the entity
         */
        inline bool
        is_external_entity( moris::moris_index aEntityIndex,
                enum EntityRank                aEntityRank ) const
        {
            if ( mFirstExtEntityInds( (moris::moris_index)aEntityRank ) <= aEntityIndex )
            {
                if ( aEntityRank == EntityRank::NODE )
                {
                    moris::moris_index tOffset         = mFirstExtEntityInds( (moris::size_t)aEntityRank );
                    moris::moris_index tNumExtEntities = mExternalEntities( (moris::size_t)aEntityRank ).size();
                    MORIS_ERROR( aEntityIndex - tOffset <= tNumExtEntities, "Requested Entity Index is out of bounds" );
                }
                return true;
            }
            else
            {
                return false;
            }
        }

        // ----------------------------------------------------------------------------------
        /*
         * Gets the external entity index
         * @param[ in ] aEntityIndex index of the entity
         * @param[ in ] aEntityRank rank of the entity
         */
        inline moris::moris_index
        get_external_entity_index(
                moris::moris_index aEntityIndex,
                enum EntityRank    aEntityRank ) const
        {
            moris::moris_index tOffset = mFirstExtEntityInds( (moris::size_t)aEntityRank );

            moris::moris_index tNumExtEntities = mExternalEntities( (moris::size_t)aEntityRank ).size();

            MORIS_ERROR( aEntityIndex - tOffset <= tNumExtEntities, "Requested Entity Index is out of bounds" );

            return aEntityIndex - mFirstExtEntityInds( (moris::moris_index)aEntityRank );
        }

        // ----------------------------------------------------------------------------------
        /**
         * Gets the owning processor of the external vertex
         * @param[ in ] aEntityIndex index of the entity
         */
        moris::moris_index
        get_external_vertex_owner( moris::moris_index aEntityIndex ) const
        {
            moris::moris_index tExternalIndex = this->get_external_entity_index( aEntityIndex, EntityRank::NODE );

            return mExternalEntities( 0 )( tExternalIndex ).get_entity_owner();
        }

        // ----------------------------------------------------------------------------------
        /**
         * A map that from entity id to local index
         * @param[ in ] aEntityIndex index of the entity
         * @param[ in ] aEntityRank rank of the entity
         */
        moris::moris_id
        get_glb_entity_id_from_entity_loc_index_external_data(
                moris::moris_id aEntityIndex,
                enum EntityRank aEntityRank ) const;

        // ----------------------------------------------------------------------------------
        /**
         * Gets nodal coordinates of nodes providing their indices
         * @param[ in ] aEntityIndex index of the entity
         */
        moris::Matrix< moris::DDRMat > const &
        get_selected_node_coordinates_loc_inds_external_data( moris::moris_index aEntityIndex ) const;

        // ----------------------------------------------------------------------------------

        void
        get_all_node_coordinates_loc_inds_external_data( moris::moris_index aStartingIndex,
                moris::Matrix< moris::DDRMat >&                             aCoordinates ) const;

        // ----------------------------------------------------------------------------------
        // TODO: [MPI] Fill in gather functions
        /**
         * Allocate entity ids
         * @param[ in ] aNumIdstoAllocate number of ids to be allocated
         * @param[ in ] aEntityRank rank of the entity
         */
        moris::moris_id
        allocate_entity_ids_external_entity_data(
                moris::size_t   aNumIdstoAllocate,
                enum EntityRank aEntityRank ) const;

        // ----------------------------------------------------------------------------------
        /**
         * Allocate entity ids
         * @param[ in ] aNumIdstoAllocate number of ids to be allocated
         * @param[ in ] aEntityRank rank of the entity
         */
        moris::Matrix< moris::IndexMat > const &
        get_local_to_global_node_map() const
        {
            return mLocalToGlobalExtNodes;
        }

        // ----------------------------------------------------------------------------------

        //                size_t
        //                capacity()
        //                {
        //                    size_t tTotal = 0;
        //                    tTotal += mFirstExtEntityInds.capacity();
        //                    tTotal += mFirstAvailableIds.capacity();
        //                    tTotal += mLocalToGlobalExtNodes.capacity();
        //                    tTotal += mFirstAvailableInds.capacity();
        //                    //tTotal += moris::internal_capacity_nested(mExternalEntities);
        //
        //                    return tTotal;
        //                }

        // ----------------------------------------------------------------------------------
    };
}    // namespace moris::mtk

#endif /* PROJECTS_MTK_SRC_CL_MTK_INTERSEC_MESH_DATA_HPP_ */
