/*
 * cl_Mesh_Parallel_Data_STK.hpp
 *
 *  Created on: Jun 20, 2017
 *      Author: ktdoble
 */

#ifndef SRC_MESH_CL_MESH_PARALLEL_DATA_STK_HPP_
#define SRC_MESH_CL_MESH_PARALLEL_DATA_STK_HPP_

// STD Includes
#include <memory>
#include <mpi.h>

// XTKL: Linear Algebra Includes
#include "mesh/cl_Mesh_Enums.hpp"

// STKL: Sierra Toolkit Mesh Includes
#include <stk_mesh/base/BulkData.hpp>    // for BulkData
#include <stk_mesh/base/MetaData.hpp>    // for MetaData
#include <stk_mesh/base/GetEntities.hpp> // for count_entities

#include "linalg/cl_XTK_Matrix.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"


using namespace xtk;

namespace mesh
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Mesh_Parallel_Data_Stk
{
public:
    /**
     * The initialization of the parallel data member variables is very long. Just think of it as a nested call to the std::vector constructor until you get to the just calling the
     * matrix constructor
     */
    Mesh_Parallel_Data_Stk(Integer const aParallelPoolSize) :
        mEntityLocaltoGlobalMap((Integer) EntityRank::END_ENUM, moris::Matrix< Integer_Matrix >(1, 1, (Integer) 0)),
        mEntitySendList((Integer) EntityRank::END_ENUM, Cell<moris::Matrix< Integer_Matrix >>(aParallelPoolSize,moris::Matrix< Integer_Matrix >(1,1,(Integer)0))),
        mEntityReceiveList((Integer)EntityRank::END_ENUM,Cell<moris::Matrix< Integer_Matrix >>(aParallelPoolSize,moris::Matrix< Integer_Matrix >(1,1,(Integer)0)))
{

}

    ~Mesh_Parallel_Data_Stk()
    {}

    void
    setup_mesh_parallel_data(stk::mesh::BulkData & aBulkData,
                             stk::mesh::MetaData & aMetaData)
    {
        setup_communication_tables_and_local_to_global_map(EntityRank::NODE, aBulkData, aMetaData);
        setup_communication_tables_and_local_to_global_map(EntityRank::EDGE, aBulkData, aMetaData);
        setup_communication_tables_and_local_to_global_map(EntityRank::FACE, aBulkData, aMetaData);
        setup_communication_tables_and_local_to_global_map(EntityRank::ELEMENT, aBulkData, aMetaData);
    }

    void update_mesh_parallel_data(Cell<Pending_Node<Real, Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes)
    {
        this->add_pending_nodes_to_local_to_global_map(aPendingNodes);
    }

    Integer get_entity_glb_id_from_entity_index(Integer & aEntityIndex,
                                                enum EntityRank aEntityRank) const
    {
        return mEntityLocaltoGlobalMap((Integer)aEntityRank)(0,aEntityIndex);
    }

    Integer get_num_of_entities_shared_with_processor_parallel_data(Integer aProcessorRank, enum EntityRank aEntityRank, bool aSendFlag) const
    {
        Integer tNumEntitiesInCommList = 0;
        if(aSendFlag)
        {
            tNumEntitiesInCommList = mEntitySendList((Integer)aEntityRank)(aProcessorRank).n_cols();
        }
        else
        {
            tNumEntitiesInCommList = mEntityReceiveList((Integer)aEntityRank)(aProcessorRank).n_cols();
        }

        return tNumEntitiesInCommList;
    }

    moris::Matrix< Integer_Matrix > const & get_local_to_global_map_parallel_data(enum EntityRank aEntityRank) const
        {
        XTK_ASSERT(aEntityRank==EntityRank::NODE,"Only allowed for nodes");
        return mEntityLocaltoGlobalMap((Integer)aEntityRank);
        }

    // Private member variables
private:
    Cell<moris::Matrix< Integer_Matrix >> mEntityLocaltoGlobalMap;
    Cell<Cell<moris::Matrix< Integer_Matrix >>> mEntitySendList;
    Cell<Cell<moris::Matrix< Integer_Matrix >>> mEntityReceiveList;

    // Private member functions
private:
    void
    setup_communication_tables_and_local_to_global_map(enum EntityRank aEntityRank,
                                                       stk::mesh::BulkData & aBulkData,
                                                       stk::mesh::MetaData const & aMetaData)
    {
        const Integer tParallelSize = aBulkData.parallel_size();
        const int tParallelRank = aBulkData.parallel_rank();

        // Declare vector of entity counts
        std::vector<uint> tEntityCounts;

        // Get all entities from meta data
        stk::mesh::Selector tSharedSelector = aMetaData.universal_part();

        // Count entities
        stk::mesh::count_entities( tSharedSelector, aBulkData, tEntityCounts );

        uint tNumEntities = static_cast<uint>(tEntityCounts[ (Integer)aEntityRank ]);

        // Resize comm lists to maximum possible
        // TODO: do this without a loop
        for(Integer i = 0; i< tParallelSize; i++)
        {
            moris::Matrix< Integer_Matrix > tSendMat(1,tNumEntities,(Integer)0);
            moris::Matrix< Integer_Matrix > tRecvMat(1,tNumEntities,(Integer)0);

            mEntitySendList((Integer)aEntityRank)(i) = tSendMat;
            mEntityReceiveList((Integer)aEntityRank)(i) = tRecvMat;
        }

        moris::Matrix< Integer_Matrix > tMapMat(1,tNumEntities,(Integer)0);
        mEntityLocaltoGlobalMap((Integer)aEntityRank)= tMapMat;

        stk::mesh::BucketVector const& shared_node_buckets =
                aBulkData.get_buckets( get_stk_entity_rank(aEntityRank) , tSharedSelector);

        Integer tCurrentIndex = 0;
        // Initialize proc counter
        // Cell #  = Proc rank
        Cell<Integer> tSendProcCounter(tParallelSize);
        Cell<Integer> tRecvProcCounter(tParallelSize);

        // Loop over shared nodes
        for(Integer i = 0; i<shared_node_buckets.size();i++)
        {
            stk::mesh::Bucket& bucket = *shared_node_buckets[i];

            //        std::cout<<bucket.size()<<std::endl;

            for(Integer j = 0; j<bucket.size(); j++)
            {
                Integer tEntityId = (Integer) aBulkData.identifier(bucket[j]);
                int tOwnerProcRank = aBulkData.parallel_owner_rank(bucket[j]);

                // Set local to global map in mesh and STK
                mEntityLocaltoGlobalMap((Integer)aEntityRank)(0, tCurrentIndex) = tEntityId;

                aBulkData.set_local_id( bucket[j], tCurrentIndex );

                std::vector<int> sharedProcs;

                // Get shared procs Ids
                aBulkData.comm_procs(aBulkData.entity_key(bucket[j]),sharedProcs);

                if(sharedProcs.size() != 0)
                {

                    // Sort (if current proc owns entity then add to send comm lists)
                    if(tOwnerProcRank == tParallelRank)
                    {
                        // loop over processors
                        for(Integer p = 0; p<sharedProcs.size(); p++)
                        {
                            if(sharedProcs[p]!=aBulkData.parallel_rank())
                            {
                                Integer tSharedProcRank = sharedProcs[p];
                                Integer tSendCount = tSendProcCounter(tSharedProcRank);
                                mEntitySendList((Integer)aEntityRank)(tSharedProcRank)(0,tSendCount) = tCurrentIndex;
                                tSendProcCounter(tSharedProcRank)++;
                            }
                        }
                    }
                    // (if current proc does not own entity then add to recv comm lists)
                    else if (tOwnerProcRank != tParallelRank)
                    {
                        for(Integer p = 0; p<sharedProcs.size(); p++)
                        {
                            if(sharedProcs[p]!=tParallelRank)
                            {
                                Integer tSharedProc = sharedProcs[p];
                                Integer tRecvCount = tRecvProcCounter(tSharedProc);
                                mEntityReceiveList((Integer)aEntityRank)(tSharedProc)(0,tRecvCount) = tCurrentIndex;
                                tRecvProcCounter(tSharedProc)++;
                            }
                        }
                    }
                }

                tCurrentIndex++;
            }
        }

        for(Integer pr = 0; pr<tParallelSize;pr++)
        {
            Integer tRecvCount = tRecvProcCounter(pr);
            Integer tSendCount = tSendProcCounter(pr);
            mEntitySendList((Integer)aEntityRank)(pr).resize(1,tSendCount);
            mEntityReceiveList((Integer)aEntityRank)(pr).resize(1,tRecvCount);
        }
    }

    stk::mesh::EntityRank get_stk_entity_rank(enum EntityRank aXtkEntityRank) const
    {
        if(aXtkEntityRank == EntityRank::NODE)
        {
            return stk::topology::NODE_RANK;
        }
        else if(aXtkEntityRank == EntityRank::EDGE)
        {
            return stk::topology::EDGE_RANK;
        }

        else if(aXtkEntityRank == EntityRank::FACE)
        {
            return stk::topology::FACE_RANK;
        }
        else if(aXtkEntityRank == EntityRank::ELEMENT)
        {
            return stk::topology::ELEMENT_RANK;
        }
        else
        {
            return stk::topology::INVALID_RANK;
        }
    }

    void add_pending_nodes_to_local_to_global_map(Cell<Pending_Node<Real, Integer, Real_Matrix, Integer_Matrix>> const & aPendingNodes)
    {
        Integer tNumNodesToAdd = aPendingNodes.size();
        Integer tPreviousSize = mEntityLocaltoGlobalMap((Integer)EntityRank::NODE).n_cols();
        Integer tNewSize = tNumNodesToAdd+tPreviousSize;
        mEntityLocaltoGlobalMap((Integer)EntityRank::NODE).resize(1,tNewSize);

        Integer tNewNodeIndex;
        Integer tNewNodeId;
        for(Integer i = 0; i< tNumNodesToAdd; i++)
        {
            tNewNodeIndex = aPendingNodes(i).get_node_index();
            tNewNodeId = aPendingNodes(i).get_node_id();
            mEntityLocaltoGlobalMap((Integer)EntityRank::NODE)(0,tNewNodeIndex) = tNewNodeId;
        }
    }
};

}

#endif /* SRC_MESH_CL_MESH_PARALLEL_DATA_STK_HPP_ */
