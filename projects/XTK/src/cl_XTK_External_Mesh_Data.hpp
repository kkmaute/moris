/*
 * cl_XTK_External_Mesh_Data.hpp
 *
 *  Created on: Jan 2, 2018
 *      Author: doble
 */

#ifndef SRC_XTK_CL_XTK_EXTERNAL_MESH_DATA_HPP_
#define SRC_XTK_CL_XTK_EXTERNAL_MESH_DATA_HPP_

#include <mpi.h>

// TPL: STK Includes

#include "cl_Matrix.hpp"
// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
#include "cl_Mesh_Enums.hpp"
// XTKL: Containers
#include "cl_Cell.hpp"

// XTKL: XTK Includes
#include "cl_XTK_Pending_Node.hpp"

#include "fn_TOL_Capacities.hpp"

// XTKL: Linear Algebra Includes
#include "cl_MPI_Tools.hpp"

using namespace moris;


#include "assert.hpp"

namespace mesh
{
    class Entity
    {
        public:

            // ----------------------------------------------------------------------------------

            Entity() :
                mGlbId(MORIS_ID_MAX),
                mLocInd(MORIS_ID_MAX),
                mOwningProc(MORIS_ID_MAX)
            {

            }

            // ----------------------------------------------------------------------------------

            ~Entity()
            {

            }

            // ----------------------------------------------------------------------------------

            void set_entity_identifiers(
                    moris::moris_id                     aGlbId,
                    moris::moris_id                     aLocInd,
                    moris::moris_id                     aOwnerProc,
                    enum moris::EntityRank              aEntityRank)
            {
                mGlbId        = aGlbId;
                mLocInd       = aLocInd;
                mOwningProc   = aOwnerProc;
                mEntityRank   = aEntityRank;
            }

            // ----------------------------------------------------------------------------------

            void set_entity_coords(moris::Matrix< moris::DDRMat > const & aCoordinates)
            {
                if (mEntityRank == moris::EntityRank::NODE)
                {
                    mEntityCoordinates = aCoordinates.copy();
                }
                else
                {
                    MORIS_ERROR(false,
                            "Entity::set_entity_coords - Only nodes should have coordinates in this context to avoid duplicate coordinate storage");
                }
            }

            // ----------------------------------------------------------------------------------

            void set_field_data(moris::Matrix< moris::DDRMat > const & aFieldData)
            {
                mNumFields = aFieldData.n_cols();
                mFieldData = aFieldData.copy();
            }

            // ----------------------------------------------------------------------------------

            moris::moris_index get_entity_loc_index() const
            {
                MORIS_ASSERT(mLocInd!=std::numeric_limits<moris::moris_id>::max(),"Index has not been set");

                return mLocInd;
            }

            // ----------------------------------------------------------------------------------

            moris::moris_index get_entity_glb_id() const
            {
                MORIS_ASSERT(mGlbId!=std::numeric_limits<moris::moris_id>::max(),"Id has not been set");
                return mGlbId;
            }

            // ----------------------------------------------------------------------------------

            moris::moris_index
            get_entity_owner() const
            {
                MORIS_ASSERT(mOwningProc!=std::numeric_limits<moris::moris_id>::max(),"Owner has not been set");
                return mOwningProc;
            }

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::DDRMat > const &
            get_entity_coords() const
            {
                return mEntityCoordinates;
            }

            // ----------------------------------------------------------------------------------

            moris::real
            get_field_data(moris::moris_index aFieldIndex) const
            {
                MORIS_ASSERT(mNumFields!=0,"Fields have not been set");
                MORIS_ASSERT(aFieldIndex<(moris::moris_index)mNumFields,"Field index is outside of bounds. Note this function should not be used directly but via STK_Mesh_Data only");
                return mFieldData(0,aFieldIndex);
            }

            // ----------------------------------------------------------------------------------

            size_t
            capacity()
            {
                size_t tTotal = 0;
                tTotal += sizeof(mGlbId);
                tTotal += sizeof(mLocInd);
                tTotal += sizeof(mOwningProc);
                tTotal += sizeof(mNumFields);
                tTotal += sizeof(mEntityRank);
                tTotal += mFieldData.capacity();
                tTotal += mEntityCoordinates.capacity();
                return tTotal;
            }

            // ----------------------------------------------------------------------------------

        private:
            moris::moris_id              mGlbId;
            moris::moris_id              mLocInd;
            moris::moris_id              mOwningProc;
            moris::size_t                mNumFields;
            enum moris::EntityRank       mEntityRank;
            moris::Matrix< moris::DDRMat > mFieldData;
            moris::Matrix< moris::DDRMat > mEntityCoordinates; // If its a node
    };
}

// ----------------------------------------------------------------------------------

namespace xtk
{
    class Mesh_External_Entity_Data
    {
        public:

            // ----------------------------------------------------------------------------------

            Mesh_External_Entity_Data() :
                mFirstExtEntityInds((moris::moris_index) EntityRank::END_ENUM, std::numeric_limits<moris::moris_index>::max()),
                mLocalToGlobalExtNodes(1,1),
                mFirstAvailableInds((moris::moris_index) EntityRank::END_ENUM, std::numeric_limits<moris::moris_index>::max()),
                mExternalEntities((moris::moris_index) EntityRank::END_ENUM, moris::Cell<mesh::Entity>(0))
            {
            }

            // ----------------------------------------------------------------------------------

            void set_up_external_entity_data(moris::mtk::Interpolation_Mesh* aMeshData)
            {
                mFirstAvailableIds.resize(7,  std::numeric_limits<moris::moris_index>::max());
                mFirstExtEntityInds.resize(7, std::numeric_limits<moris::moris_index>::max());
                mFirstAvailableInds.resize(7, std::numeric_limits<moris::moris_index>::max());

                int tProcessorRank = moris::par_rank();

                moris::moris_id tFirstNode = aMeshData->get_max_entity_id(EntityRank::NODE)+1;
                //        moris::moris_id tFirstEdge = aMeshData->get_max_entity_id(EntityRank::EDGE);
                //        moris::moris_id tFirstFace = aMeshData->get_max_entity_id(EntityRank::FACE);
                moris::moris_id tFirstElem = aMeshData->get_max_entity_id(EntityRank::ELEMENT)+1;
                
                if(tProcessorRank == 0)
                {
                    // Processor 1 (rank 0) is responsible for first available Ids
                    mFirstAvailableIds(0) = tFirstNode;
                    mFirstAvailableIds(1) = MORIS_ID_MAX;
                    mFirstAvailableIds(2) = MORIS_ID_MAX;
                    mFirstAvailableIds(3) = tFirstElem;
                }

                mFirstExtEntityInds(0) = aMeshData->get_num_entities(moris::EntityRank::NODE);
                mFirstExtEntityInds(1) = aMeshData->get_num_entities(moris::EntityRank::EDGE);
                mFirstExtEntityInds(2) = aMeshData->get_num_entities(moris::EntityRank::FACE);
                mFirstExtEntityInds(3) = aMeshData->get_num_entities(moris::EntityRank::ELEMENT);

                mFirstAvailableInds(0) = mFirstExtEntityInds(0);
                mFirstAvailableInds(1) = mFirstExtEntityInds(1);
                mFirstAvailableInds(2) = mFirstExtEntityInds(2);
                mFirstAvailableInds(3) = mFirstExtEntityInds(3);
            }

            // ----------------------------------------------------------------------------------

            moris::moris_index
            get_first_available_index_external_data(enum EntityRank aEntityRank) const
            {
                return mFirstAvailableInds((moris::moris_index)aEntityRank);
            }

            // ----------------------------------------------------------------------------------

            void update_first_available_index_external_data(moris::moris_index aNewFirstAvailableIndex, enum EntityRank aEntityRank)
            {
                mFirstAvailableInds((moris::moris_index)aEntityRank) = aNewFirstAvailableIndex;
            }

            // ----------------------------------------------------------------------------------

            void
            batch_create_new_nodes_external_data(
                    Cell<moris_index>                    const & aNewNodeIds,
                    Cell<moris_index>                    const & aNewNodeIndices,
                    Cell<moris_index>                    const & aNewNodeOwners,
                    Cell<moris::Matrix< moris::DDRMat >> const & aNewNodeCoordinates)
            {
                moris::moris_index tEntRankInd  = (moris::moris_index)EntityRank::NODE;
                moris::size_t      tAddSize     = aNewNodeIds.size();
                moris::size_t      tInitialSize = mExternalEntities(tEntRankInd).size();

                // Initialize
                moris::size_t j    = 0;
                moris::moris_index tInd   = MORIS_INDEX_MAX;
                moris::moris_id    tId    = MORIS_ID_MAX;
                moris::moris_id    tOwner = MORIS_ID_MAX;

                // Resize
                mExternalEntities(tEntRankInd).resize(tInitialSize+tAddSize,mesh::Entity());
                mLocalToGlobalExtNodes.resize(1,tInitialSize+tAddSize);

                for(moris::size_t i = tInitialSize; i<tAddSize+tInitialSize;i++)
                {
                    // Add information to entities
                    tInd    = aNewNodeIndices(j);
                    tId     = aNewNodeIds(j);
                    tOwner  = aNewNodeOwners(j);

                    mLocalToGlobalExtNodes(tInd-mFirstExtEntityInds(0)) = tId;
                    mExternalEntities(tEntRankInd)(tInd-mFirstExtEntityInds(0)).set_entity_identifiers(tId,tInd,tOwner,moris::EntityRank::NODE);
                    mExternalEntities(tEntRankInd)(tInd-mFirstExtEntityInds(0)).set_entity_coords(aNewNodeCoordinates(j));
                    j++;
                }
            }

            // ----------------------------------------------------------------------------------

            void
            batch_create_new_nodes_external_data(
                    moris::Matrix< moris::IndexMat >    const & aNewNodeIds,
                    moris::Matrix< moris::IndexMat >    const & aNewNodeIndices,
                    moris::Matrix< moris::IndexMat >    const & aNewNodeOwners,
                    moris::Matrix< moris::DDRMat >      const & aNewNodeCoordinates)
            {
                moris::moris_index tEntRankInd  = (moris::moris_index)EntityRank::NODE;
                moris::size_t      tAddSize     = aNewNodeIds.numel();
                moris::size_t      tInitialSize = mExternalEntities(tEntRankInd).size();

                // Initialize
                moris::size_t j    = 0;
                moris::moris_index tInd   = MORIS_INDEX_MAX;
                moris::moris_id    tId    = MORIS_ID_MAX;
                moris::moris_id    tOwner = MORIS_ID_MAX;

                // Resize
                mExternalEntities(tEntRankInd).resize((tInitialSize+tAddSize),mesh::Entity());
                mLocalToGlobalExtNodes.resize(1,(tInitialSize+tAddSize));

                for(moris::size_t i = tInitialSize; i<tAddSize+tInitialSize;i++)
                {
                    // Add information to entities
                    tInd    = aNewNodeIndices(j);
                    tId     = aNewNodeIds(j);
                    tOwner  = aNewNodeOwners(j);
                    mLocalToGlobalExtNodes(tInd-mFirstExtEntityInds(0)) = tId;
                    mExternalEntities(tEntRankInd)(i).set_entity_identifiers(tId,tInd,tOwner,moris::EntityRank::NODE);
                    mExternalEntities(tEntRankInd)(i).set_entity_coords(aNewNodeCoordinates.get_row(j));
                    j++;
                }
            }

            // ----------------------------------------------------------------------------------

            moris::size_t get_num_entities_external_data(enum EntityRank aEntityRank) const
            {
                return mExternalEntities((moris::size_t)aEntityRank).size();
            }

            // ----------------------------------------------------------------------------------

            inline
            bool is_external_entity(moris::moris_index aEntityIndex,
                    enum EntityRank aEntityRank) const
            {
                if(mFirstExtEntityInds((moris::moris_index)aEntityRank)<=aEntityIndex)
                {
                    if(aEntityRank == EntityRank::NODE)
                    {
                        moris::moris_index tOffset = mFirstExtEntityInds((moris::size_t)aEntityRank);
                        moris::moris_index tNumExtEntities = mExternalEntities((moris::size_t)aEntityRank).size();
                        MORIS_ERROR(aEntityIndex-tOffset<=tNumExtEntities,"Requested Entity Index is out of bounds");
                    }
                    return true;
                }
                else
                {
                    return false;
                }
            }

            // ----------------------------------------------------------------------------------

            inline
            moris::moris_index
            get_external_entity_index(
                    moris::moris_index aEntityIndex,
                    enum EntityRank    aEntityRank) const
            {
                moris::moris_index tOffset = mFirstExtEntityInds((moris::size_t)aEntityRank);

                moris::moris_index tNumExtEntities = mExternalEntities((moris::size_t)aEntityRank).size();

                MORIS_ERROR(aEntityIndex-tOffset<=tNumExtEntities,"Requested Entity Index is out of bounds");

                return  aEntityIndex - mFirstExtEntityInds((moris::moris_index)aEntityRank);
            }

            // ----------------------------------------------------------------------------------

            moris::moris_index
            get_external_vertex_owner(moris::moris_index aEntityIndex) const
            {
                moris::moris_index tExternalIndex = this->get_external_entity_index(aEntityIndex,EntityRank::NODE);

                return mExternalEntities(0)(tExternalIndex).get_entity_owner();
            }

            // ----------------------------------------------------------------------------------

            moris::moris_id get_glb_entity_id_from_entity_loc_index_external_data(moris::moris_id aEntityIndex, enum EntityRank aEntityRank) const
            {
                moris::size_t tEntityRankIndex = (moris::size_t)aEntityRank;

                moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

                return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_glb_id();
            }

            // ----------------------------------------------------------------------------------

            moris::Matrix< moris::DDRMat > const &
            get_selected_node_coordinates_loc_inds_external_data(moris::moris_index aEntityIndex) const
            {
                moris::size_t tEntityRankIndex = (moris::size_t)EntityRank::NODE;
                moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

                return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_coords();
            }

            // ----------------------------------------------------------------------------------

            void
            get_all_node_coordinates_loc_inds_external_data(moris::moris_index aStartingIndex,
                    moris::Matrix< moris::DDRMat > & aCoordinates) const
            {
                moris::size_t tEntityRankIndex = (moris::size_t)EntityRank::NODE;
                moris::size_t tNumNodes = this->get_num_entities_external_data(EntityRank::NODE);

                for( moris::moris_index i = 0; i<(moris::moris_index)tNumNodes; i++)
                {
                    const moris::Matrix< moris::DDRMat > & tCoordinateRow = mExternalEntities(tEntityRankIndex)(i).get_entity_coords();
                    aCoordinates.set_row(aStartingIndex,tCoordinateRow);
                    aStartingIndex++;
                }
            }

            // ----------------------------------------------------------------------------------

            //TODO: [MPI] Fill in gather functions
            moris::moris_id
            allocate_entity_ids_external_entity_data(
                    moris::size_t aNumIdstoAllocate,
                    enum EntityRank aEntityRank) const
            {
                int tProcRank = par_rank();
                int tProcSize = par_size();

                // size_t is defined as uint here because of aNumRequested
                //Initialize gathered information outputs (information which will be scattered across processors)
                moris::Cell<moris::moris_id> aGatheredInfo;
                moris::Cell<moris::moris_id> tFirstId(1);
                moris::Cell<moris::moris_id> tNumIdsRequested(1);

                tNumIdsRequested(0) = (moris::moris_id)aNumIdstoAllocate;

                moris::gather(tNumIdsRequested,aGatheredInfo);

                moris::Cell<moris::moris_id> tProcFirstID(tProcSize);

                moris::moris_id tEntityRankIndex = (moris::moris_id)aEntityRank;

                if(tProcRank == 0)
                {
                    // Loop over entities print the number of entities requested by each processor
                    for (int iProc = 0; iProc < tProcSize; ++iProc)
                    {
                        // Give each processor their desired amount of IDs
                        tProcFirstID(iProc) = mFirstAvailableIds(tEntityRankIndex);

                        // Increment the first available node ID
                        mFirstAvailableIds(tEntityRankIndex) = mFirstAvailableIds(tEntityRankIndex)+aGatheredInfo(iProc);
                    }
                }

                moris::scatter(tProcFirstID,tFirstId);

                return tFirstId(0);
            }

            // ----------------------------------------------------------------------------------

            moris::real get_entity_field_value_external_data(
                    moris::moris_index  aEntityIndex,
                    enum EntityRank     aFieldEntityRank,
                    std::string const & aFieldName) const
            {
                moris::size_t tFieldIndex = this->get_field_index(aFieldName);
                moris::size_t tEntityRankIndex = (moris::size_t)aFieldEntityRank;
                moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);
                return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_field_data(tFieldIndex);
            }

            // ----------------------------------------------------------------------------------

            moris::Matrix<moris::IndexMat> const &
            get_local_to_global_node_map() const
            {
                return mLocalToGlobalExtNodes;
            }

            // ----------------------------------------------------------------------------------

            size_t
            capacity()
            {
                size_t tTotal = 0;
                tTotal += mFirstExtEntityInds.capacity();
                tTotal += mFirstAvailableIds.capacity();
                tTotal += mLocalToGlobalExtNodes.capacity();
                tTotal += mFirstAvailableInds.capacity();
                tTotal += moris::internal_capacity_nested(mExternalEntities);
                tTotal += moris::internal_capacity(mFieldNames);

                return tTotal;
            }

            // ----------------------------------------------------------------------------------

        private:

            moris::Cell<moris::moris_index> mFirstExtEntityInds;

            // Owned by proc rank 0, other procs UINT_MAX
            // Mutable to preserve const in the allocate entity ids function
            mutable moris::Cell<moris::moris_id> mFirstAvailableIds;

            // Local to Global Node Map
            moris::Matrix<moris::IdMat> mLocalToGlobalExtNodes;

            // Each processor tracks this value
            moris::Cell<moris::moris_index> mFirstAvailableInds;

            // Entity Rank outside, then entity objects inside
            moris::Cell<moris::Cell<mesh::Entity>>mExternalEntities;

            // Fields
            moris::Cell<std::string> mFieldNames;

        private:

            // ----------------------------------------------------------------------------------

            void register_fields(moris::Cell<std::string> const & aFieldNames)
            {
                if(mFieldNames.size()!=0)
                {
                    // Check to see if they match
                    for(moris::size_t i = 0 ;i<mFieldNames.size(); i++)
                    {
                        if(mFieldNames(i).compare(aFieldNames(i))!=0)
                        {
                            bool tBreak = true;
                            MORIS_ERROR(!tBreak,"Fields provided do not match the ones already in external mesh data. Currently cannot add fields to external data");
                        }
                    }
                }
                else
                {
                    mFieldNames = aFieldNames;
                }
            }

            // ----------------------------------------------------------------------------------

            moris::size_t get_field_index(std::string const & aFieldName) const
            {
                bool tSuccess = false;
                moris::size_t tFieldIndex = 0;
                for(moris::size_t i = 0 ;i<mFieldNames.size(); i++)
                {
                    if(mFieldNames(i).compare(aFieldName) == 0)
                    {
                        tSuccess = true;
                        tFieldIndex = i;
                        break;
                    }
                }

                MORIS_ERROR(tSuccess,"Could not locate the specified field name");
                return tFieldIndex;
            }

            // ----------------------------------------------------------------------------------
    };
}

#endif /* SRC_XTK_CL_XTK_EXTERNAL_MESH_DATA_HPP_ */
