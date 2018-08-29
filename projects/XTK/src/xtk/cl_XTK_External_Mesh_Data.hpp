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
#include <stk_mesh/base/BulkData.hpp>    // for BulkData

#include "linalg/cl_XTK_Matrix_Base.hpp"
// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Enums.hpp"
#include "mesh/cl_Mesh_Entity.hpp"
// XTKL: Containers
#include "containers/cl_XTK_Cell.hpp"

// XTKL: XTK Includes
#include "xtk/cl_XTK_Pending_Node.hpp"

// XTKL: Linear Algebra Includes
#include "tools/cl_MPI_Tools.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Mesh_External_Entity_Data
{
public:
    Integer INTEGER_MAX = std::numeric_limits<Integer>::max();
    Mesh_External_Entity_Data() :
            mFirstExtEntityInds((Integer) EntityRank::END_ENUM, INTEGER_MAX),
            mFirstAvailableIds((Integer) EntityRank::END_ENUM, INTEGER_MAX),
            mFirstAvailableInds((Integer) EntityRank::END_ENUM, INTEGER_MAX),
            mExternalEntities((Integer) EntityRank::END_ENUM, xtk::Cell<mesh::Entity<Real, Integer, Real_Matrix>>(0))
    {
    }

    void set_up_external_entity_data(Integer aNumNodes,
                                     Integer aNumEdges,
                                     Integer aNumFaces,
                                     Integer aNumElements,
                                     stk::mesh::BulkData & aBulkData)
    {
        mFirstAvailableIds.resize(4,  INTEGER_MAX);
        mFirstExtEntityInds.resize(4, INTEGER_MAX);
        mFirstAvailableInds.resize(4, INTEGER_MAX);

        int tProcessorRank;
        MPI_Comm_rank(MPI_COMM_WORLD, &tProcessorRank);

        std::vector<stk::mesh::EntityId> tFirstNode;
        std::vector<stk::mesh::EntityId> tFirstEdge;
        std::vector<stk::mesh::EntityId> tFirstFace;
        std::vector<stk::mesh::EntityId> tFirstElem;
        aBulkData.generate_new_ids(stk::topology::NODE_RANK, 1, tFirstNode);
        aBulkData.generate_new_ids(stk::topology::EDGE_RANK, 1, tFirstEdge);
        aBulkData.generate_new_ids(stk::topology::FACE_RANK, 1, tFirstFace);
        aBulkData.generate_new_ids(stk::topology::ELEMENT_RANK, 1, tFirstElem);


        if(tProcessorRank == 0)
        {
            // Processor 1 (rank 0) is responsible for first available Ids
            mFirstAvailableIds(0) = (Integer) tFirstNode[0];
            mFirstAvailableIds(1) = (Integer) tFirstEdge[0];
            mFirstAvailableIds(2) = (Integer) tFirstFace[0];
            mFirstAvailableIds(3) = (Integer) tFirstElem[0];
        }

        mFirstExtEntityInds(0) = aNumNodes;
        mFirstExtEntityInds(1) = aNumEdges;
        mFirstExtEntityInds(2) = aNumFaces;
        mFirstExtEntityInds(3) = aNumElements;

        mFirstAvailableInds(0) = mFirstExtEntityInds(0);
        mFirstAvailableInds(1) = mFirstExtEntityInds(1);
        mFirstAvailableInds(2) = mFirstExtEntityInds(2);
        mFirstAvailableInds(3) = mFirstExtEntityInds(3);
    }

    Integer
    get_first_available_index_external_data(enum EntityRank aEntityRank) const
    {
        return mFirstAvailableInds((Integer)aEntityRank);
    }

    void update_first_available_index_external_data(Integer aNewFirstAvailableIndex, enum EntityRank aEntityRank)
    {
        mFirstAvailableInds((Integer)aEntityRank) = aNewFirstAvailableIndex;
    }

    void batch_create_new_nodes_external_data(xtk::Cell<xtk::Pending_Node<Real,Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes,
                                              bool aFieldsToAdd=false)
    {
        Integer INTEGER_MAX = std::numeric_limits<Integer>::max();

        Integer tEntInd       = (Integer)EntityRank::NODE;
        Integer tAddSize      = aPendingNodes.size();
        Integer tInitialSize  = mExternalEntities(tEntInd).size();

        Integer j    = 0;
        Integer tInd = INTEGER_MAX;
        Integer tId  = INTEGER_MAX;

        // Resize
        mExternalEntities(tEntInd).resize((tInitialSize+tAddSize),mesh::Entity<Real,Integer,Real_Matrix>());

        for(Integer i = tInitialSize; i<tAddSize+tInitialSize;i++)
        {
            // Add information to entities
            tInd    = aPendingNodes(j).get_node_index();
            tId     = aPendingNodes(j).get_node_id();
            moris::Mat_New<Real, Real_Matrix> const & tCoords = aPendingNodes(j).get_coordinates();
            mExternalEntities(tEntInd)(i).set_entity_identifiers(tId,tInd,EntityRank::NODE);
            mExternalEntities(tEntInd)(i).set_entity_coords(tCoords);
            if(aFieldsToAdd)
            {
                moris::Mat_New<Real, Real_Matrix> const & tFieldData = aPendingNodes(j).get_field_data();
                mExternalEntities(tEntInd)(i).set_field_data(tFieldData);
            }
            j++;
        }
    }

    void batch_create_new_nodes_with_fields_external_data(xtk::Cell<xtk::Pending_Node<Real,Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes,
                                                  xtk::Cell<std::string> const & aFieldNames)
        {

        this->register_fields(aFieldNames);
        this->batch_create_new_nodes_external_data(aPendingNodes,true);
        }


    Integer get_num_entities_external_data(enum EntityRank aEntityRank) const
    {
        return mExternalEntities((Integer)aEntityRank).size();
    }

    bool is_external_entity(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        if(mFirstExtEntityInds((Integer)aEntityRank)<=aEntityIndex)
        {
            Integer tOffset = mFirstExtEntityInds((Integer)aEntityRank);
            Integer tNumExtEntities = mExternalEntities((Integer)aEntityRank).size();
            XTK_ASSERT(aEntityIndex-tOffset<=tNumExtEntities,"Requested Entity Index is out of bounds");
            return true;
        }
        else
        {
            return false;
        }
    }

    Integer get_glb_entity_id_from_entity_loc_index_external_data(Integer aEntityIndex, enum EntityRank aEntityRank) const
    {
        Integer tEntityRankIndex = (Integer)aEntityRank;
        Integer tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_glb_id();
    }


    moris::Mat_New<Real, Real_Matrix> const &
    get_selected_node_coordinates_loc_inds_external_data(Integer aEntityIndex) const
    {
        Integer tEntityRankIndex = (Integer)EntityRank::NODE;
        Integer tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);

        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_entity_coords();
    }


    void
    get_all_node_coordinates_loc_inds_external_data(Integer aStartingIndex, moris::Mat_New<Real, Real_Matrix> & aCoordinates) const
    {
        Integer tEntityRankIndex = (Integer)EntityRank::NODE;
        Integer tNumNodes = this->get_num_entities_external_data(EntityRank::NODE);

        for( Integer i = 0; i<tNumNodes; i++)
        {
            const moris::Mat_New<Real, Real_Matrix> & tCoordinateRow = mExternalEntities(tEntityRankIndex)(i).get_entity_coords();
            aCoordinates.set_row(aStartingIndex,tCoordinateRow);
            aStartingIndex++;
        }
    }

    //TODO: [MPI] Fill in gather functions
    Integer allocate_entity_ids_external_entity_data(Integer aNumIdstoAllocate, enum EntityRank aEntityRank) const
    {
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
        MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

        // size_t is defined as uint here because of aNumRequested
        //Initialize gathered information outputs (information which will be scattered across processors)
        xtk::Cell<Integer> aGatheredInfo;
        xtk::Cell<Integer> tFirstId(1);
        xtk::Cell<Integer> tNumIdsRequested(1);

        tNumIdsRequested(0) = (xtk::uint)aNumIdstoAllocate;

        xtk::gather(tNumIdsRequested,aGatheredInfo);

        xtk::Cell<Integer> tProcFirstID(tProcSize);

        Integer tEntityRankIndex = (Integer)aEntityRank;

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

    Real get_entity_field_value_external_data(Integer aEntityIndex,
                                              enum EntityRank aFieldEntityRank,
                                              std::string const & aFieldName) const
    {
        Integer tFieldIndex = this->get_field_index(aFieldName);
        Integer tEntityRankIndex = (Integer)aFieldEntityRank;
        Integer tExternalIndex = aEntityIndex - mFirstExtEntityInds(tEntityRankIndex);
        return mExternalEntities(tEntityRankIndex)(tExternalIndex).get_field_data(tFieldIndex);
    }

private:
    xtk::Cell<Integer> mFirstExtEntityInds;

    // Owned by proc rank 0, other procs UINT_MAX
    // Mutable to preserve const in the allocate entity ids function
    mutable xtk::Cell<Integer> mFirstAvailableIds;

    // Each processor tracks this value
    xtk::Cell<Integer> mFirstAvailableInds;

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
