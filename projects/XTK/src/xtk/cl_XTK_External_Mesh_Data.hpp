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

#include "linalg/cl_XTK_Matrix_Base.hpp"
// XTKL: Mesh Includes
#include "cl_MTK_Mesh.hpp"
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Entity.hpp"
// XTKL: Containers
#include "containers/cl_XTK_Cell.hpp"

// XTKL: XTK Includes
#include "xtk/cl_XTK_Pending_Node.hpp"

// XTKL: Linear Algebra Includes
#include "tools/cl_MPI_Tools.hpp"

using namespace moris;

namespace xtk
{
template<typename Real,
        typename Integer,
        typename Real_Matrix,
        typename Integer_Matrix>
class Mesh_External_Entity_Data
{
public:
    Mesh_External_Entity_Data() :
            mFirstExtEntityInds((moris::moris_index) EntityRank::END_ENUM, std::numeric_limits<moris::moris_index>::max()),
            mLocalToGlobalExtNodes(1,1),
            mFirstAvailableInds((moris::moris_index) EntityRank::END_ENUM, std::numeric_limits<moris::moris_index>::max()),
            mExternalEntities((moris::moris_index) EntityRank::END_ENUM, xtk::Cell<mesh::Entity<Real, Integer, Real_Matrix>>(0))
    {
    }

    void set_up_external_entity_data(moris::mtk::Mesh* & aMeshData)
    {
        mFirstAvailableIds.resize(4,  std::numeric_limits<moris::moris_index>::max());
        mFirstExtEntityInds.resize(4, std::numeric_limits<moris::moris_index>::max());
        mFirstAvailableInds.resize(4, std::numeric_limits<moris::moris_index>::max());

        int tProcessorRank;
        MPI_Comm_rank(get_comm(), &tProcessorRank);


        moris::moris_id tFirstNode = aMeshData->generate_unique_entity_ids(1,moris::EntityRank::NODE)(0);
        moris::moris_id tFirstEdge = aMeshData->generate_unique_entity_ids(1,moris::EntityRank::EDGE)(0);
        moris::moris_id tFirstFace = aMeshData->generate_unique_entity_ids(1,moris::EntityRank::FACE)(0);
        moris::moris_id tFirstElem = aMeshData->generate_unique_entity_ids(1,moris::EntityRank::ELEMENT)(0);

        if(tProcessorRank == 0)
        {
            // Processor 1 (rank 0) is responsible for first available Ids
            mFirstAvailableIds(0) = tFirstNode;
            mFirstAvailableIds(1) = tFirstEdge;
            mFirstAvailableIds(2) = tFirstFace;
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


    moris::moris_index
    get_first_available_index_external_data(enum EntityRank aEntityRank) const
    {
        return mFirstAvailableInds((moris::moris_index)aEntityRank);
    }

    void update_first_available_index_external_data(moris::moris_index aNewFirstAvailableIndex, enum EntityRank aEntityRank)
    {
        mFirstAvailableInds((moris::moris_index)aEntityRank) = aNewFirstAvailableIndex;
    }


    void
    batch_create_new_nodes_external_data(
            Cell<Pending_Node<Real,Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes)
    {
        moris::moris_index Index_Max = std::numeric_limits<moris::moris_index>::max();

        moris::moris_index tEntInd = (moris::moris_index)EntityRank::NODE;
        Integer tAddSize      = aPendingNodes.size();
        Integer tInitialSize  = mExternalEntities(tEntInd).size();

        Integer j    = 0;
        moris::moris_index tInd = Index_Max;
        moris::moris_id    tId  = Index_Max;

        // Resize
        mExternalEntities(tEntInd).resize((tInitialSize+tAddSize),mesh::Entity<Real,Integer,Real_Matrix>());
        mLocalToGlobalExtNodes.resize(1,(tInitialSize+tAddSize));

        for(Integer i = tInitialSize; i<tAddSize+tInitialSize;i++)
        {
            // Add information to entities
            tInd    = aPendingNodes(j).get_node_index();
            tId     = aPendingNodes(j).get_node_id();
            mLocalToGlobalExtNodes(tInd-mFirstExtEntityInds(0)) = tId;
            moris::Matrix< Real_Matrix > const & tCoords = aPendingNodes(j).get_coordinates();
            mExternalEntities(tEntInd)(i).set_entity_identifiers(tId,tInd,EntityRank::NODE);
            mExternalEntities(tEntInd)(i).set_entity_coords(tCoords);
            j++;
        }
    }



    Integer get_num_entities_external_data(enum EntityRank aEntityRank) const
    {
        return mExternalEntities((Integer)aEntityRank).size();
    }

    inline
    bool is_external_entity(moris::moris_index aEntityIndex,
                            enum EntityRank aEntityRank) const
    {
        if(mFirstExtEntityInds((moris::moris_index)aEntityRank)<=aEntityIndex)
        {
            moris::moris_index tOffset = mFirstExtEntityInds((Integer)aEntityRank);
            moris::moris_index tNumExtEntities = mExternalEntities((Integer)aEntityRank).size();
            XTK_ASSERT(aEntityIndex-tOffset<=tNumExtEntities,"Requested Entity Index is out of bounds");
            return true;
        }
        else
        {
            return false;
        }
    }

    inline
    moris::moris_index
    get_external_entity_index(moris::moris_index aEntityIndex,
                              enum EntityRank    aEntityRank) const
    {
        return  aEntityIndex - mFirstExtEntityInds((moris::moris_index)aEntityRank);
    }

    moris::moris_id get_glb_entity_id_from_entity_loc_index_external_data(moris::moris_id aEntityIndex, enum EntityRank aEntityRank) const
    {
        Integer tEntityRankIndex = (Integer)aEntityRank;
        moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_glb_id();
    }


    moris::Matrix< Real_Matrix > const &
    get_selected_node_coordinates_loc_inds_external_data(moris::moris_index aEntityIndex) const
    {
        Integer tEntityRankIndex = (Integer)EntityRank::NODE;
        moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_coords();
    }


    void
    get_all_node_coordinates_loc_inds_external_data(moris::moris_index aStartingIndex,
                                                    moris::Matrix< Real_Matrix > & aCoordinates) const
    {
        Integer tEntityRankIndex = (Integer)EntityRank::NODE;
        Integer tNumNodes = this->get_num_entities_external_data(EntityRank::NODE);

        for( moris::moris_index i = 0; i<(moris::moris_index)tNumNodes; i++)
        {
            const moris::Matrix< Real_Matrix > & tCoordinateRow = mExternalEntities(tEntityRankIndex)(i).get_entity_coords();
            aCoordinates.set_row(aStartingIndex,tCoordinateRow);
            aStartingIndex++;
        }
    }

    //TODO: [MPI] Fill in gather functions
    moris::moris_id
    allocate_entity_ids_external_entity_data(
            Integer aNumIdstoAllocate,
            enum EntityRank aEntityRank) const
    {
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
        MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

        // size_t is defined as uint here because of aNumRequested
        //Initialize gathered information outputs (information which will be scattered across processors)
        xtk::Cell<moris::moris_id> aGatheredInfo;
        xtk::Cell<moris::moris_id> tFirstId(1);
        xtk::Cell<moris::moris_id> tNumIdsRequested(1);

        tNumIdsRequested(0) = (moris::moris_id)aNumIdstoAllocate;

        xtk::gather(tNumIdsRequested,aGatheredInfo);

        xtk::Cell<moris::moris_id> tProcFirstID(tProcSize);

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

        xtk::scatter(tProcFirstID,tFirstId);


        return tFirstId(0);
    }

    Real get_entity_field_value_external_data(moris::moris_index aEntityIndex,
                                              enum EntityRank aFieldEntityRank,
                                              std::string const & aFieldName) const
    {
        Integer tFieldIndex = this->get_field_index(aFieldName);
        Integer tEntityRankIndex = (Integer)aFieldEntityRank;
        moris::moris_index tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);
        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_field_data(tFieldIndex);
    }

    moris::Matrix<moris::IndexMat> const &
    get_local_to_global_node_map() const
    {
        return mLocalToGlobalExtNodes;
    }

private:
    xtk::Cell<moris::moris_index> mFirstExtEntityInds;

    // Owned by proc rank 0, other procs UINT_MAX
    // Mutable to preserve const in the allocate entity ids function
    mutable xtk::Cell<moris::moris_id> mFirstAvailableIds;


    // Local to Global Node Map
    moris::Matrix<moris::IdMat> mLocalToGlobalExtNodes;


    // Each processor tracks this value
    xtk::Cell<moris::moris_index> mFirstAvailableInds;

    // Entity Rank outside, then entity objects inside
    xtk::Cell<xtk::Cell<mesh::Entity<Real,Integer,Real_Matrix>>>mExternalEntities;

    // Fields
    xtk::Cell<std::string> mFieldNames;



private:
    void register_fields(xtk::Cell<std::string> const & aFieldNames)
    {
        if(mFieldNames.size()!=0)
        {
            // Check to see if they match
            bool tFieldsMatch = true;
            for(Integer i = 0 ;i<mFieldNames.size(); i++)
            {
                if(mFieldNames(i).compare(aFieldNames(i))!=0)
                {
                    bool tBreak = true;
                    XTK_ASSERT(!tBreak,"Fields provided do not match the ones already in external mesh data. Currently cannot add fields to external data");
                }
            }

        }

        else
        {
            mFieldNames = aFieldNames;
        }
    }

    Integer get_field_index(std::string const & aFieldName) const
    {
        bool tSuccess = false;
        Integer tFieldIndex = 0;
        for(Integer i = 0 ;i<mFieldNames.size(); i++)
        {
            if(mFieldNames(i).compare(aFieldName) == 0)
            {
                tSuccess = true;
                tFieldIndex = i;
                break;
            }
        }

        XTK_ASSERT(tSuccess,"Could not locate the specified field name");
        return tFieldIndex;
    }

};
}

#endif /* SRC_XTK_CL_XTK_EXTERNAL_MESH_DATA_HPP_ */
