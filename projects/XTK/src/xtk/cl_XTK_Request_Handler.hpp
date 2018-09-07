/*
 * cl_XTK_Request_Handler.hpp
 *
 *  Created on: Jul 1, 2017
 *      Author: ktdoble
 */

#ifndef SRC_XTK_CL_XTK_REQUEST_HANDLER_HPP_
#define SRC_XTK_CL_XTK_REQUEST_HANDLER_HPP_

#include <limits>
#include <mpi.h>

#include "tools/fn_Pairing.hpp"

// XTKL: Linear Algebra Includes
#include "linalg/cl_XTK_Matrix.hpp"

// XTKL: Mesh Includes
#include "mesh/cl_Mesh_Data.hpp"
#include "mesh/cl_Mesh_Enums.hpp"

// XTKL: Xtk includes
#include "xtk/cl_XTK_Pending_Node.hpp"
#include "xtk/cl_XTK_Active_Process_Manager.hpp"
#include "xtk/cl_XTK_Entity_Tracker.hpp"
#include "xtk/cl_XTK_Cut_Mesh.hpp"

// XTKL: Container includes
#include "containers/cl_XTK_Cell.hpp"

// XTKL: Logging and Assertion Includes
#include "assert/fn_xtk_assert.hpp"
#include "ios/cl_Logger.hpp"

// XTKL: Tools includes
#include "tools/cl_MPI_Tools.hpp"

// Topology
#include "topology/cl_XTK_Topology.hpp"

namespace xtk
{
template<typename Real, typename Integer, typename Real_Matrix, typename Integer_Matrix>
class Request_Handler
{
    /**
     * RequestHandler stores all of the requests for creating entities of same
     * parent entity rank and child entity rank
     */
    // Forward declaration
    Integer INTEGER_MAX = std::numeric_limits<Integer>::max();

public:
    Request_Handler(Integer aNumExpectedRequests,
                    Integer aNumChildrenAllowed,
                    enum EntityRank aParentEntityRank,
                    enum EntityRank aChildEntityRank,
                    mesh::Mesh_Data<Real, Integer,Real_Matrix, Integer_Matrix> & aParentMesh,
                    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aXTKMesh) :
    mRequestCounter(0),
    mNumChildrenAllowed(aNumChildrenAllowed),
    mEntityTracker(aParentEntityRank, aChildEntityRank, aParentMesh.get_num_entities(aParentEntityRank), aNumChildrenAllowed),
    mParentEntityRank(aParentEntityRank),
    mChildEntityRank(aChildEntityRank),
    mEntityRequestInfo(aNumExpectedRequests, 2, INTEGER_MAX),
    mMesh(aParentMesh),
    mXTKMesh(aXTKMesh)
    {
        if (aChildEntityRank == EntityRank::NODE)
        {
            mPendingNodes.resize(aNumExpectedRequests, Pending_Node<Real, Integer,Real_Matrix, Integer_Matrix>());
        }
    }

    ~Request_Handler()
    {

    }


    /*
     * Takes a parent entity index and returns whether or not the provided node has been create
     */
    bool
    has_parent_entity_been_used(Integer aParentEntityIndex)
    {
        return mEntityTracker.is_parent_entity_used(aParentEntityIndex);
    }

    bool
    has_parent_entity_been_used(Integer aParentEntityIndex,
                                Integer aSecondaryEntityIdentifier)
    {
//        if(pIdInd(0) == NULL)
//        {
//            return false;
//        }


        return false;
    }

    /*
     * When you set a request, you get a pointer to the location where
     * the node index is going to be placed
     */
    Integer*
    set_request_info(Integer aParentEntityIndex,
                     Topology<Real,Integer,Real_Matrix,Integer_Matrix> const & aParentTopology,
                     moris::Matrix<Real,Real_Matrix> const &                    aChildCoordsGlb,
                     moris::Matrix<Real,Real_Matrix> const &                    aChildCoordsLoc)
    {
        // Ask entity tracker if this parent entity index has been used yet
        bool tUse = mEntityTracker.is_parent_entity_used(aParentEntityIndex);

        XTK_ASSERT(aChildCoordsGlb.n_rows() == 1, "Coordinates submitted to set_request_info need to be one row (x,y,z)");

        if (tUse != true)
        {
            if(mRequestCounter>=mEntityRequestInfo.n_rows())
            {
                XTK_ERROR<<"Warning: Not enough space allocated by constructor of Request handler, check aNumExpectedRequests. Dynamic Allocation"<<std::endl;
                Integer tNumRows  = mEntityRequestInfo.n_rows();
                Integer tNumCols  = mEntityRequestInfo.n_cols();
                mEntityRequestInfo.resize(2*tNumRows, tNumCols);
                mPendingNodes.resize(2*tNumRows, Pending_Node<Real, Integer,Real_Matrix, Integer_Matrix>());

            }

            mEntityTracker.mark_entity_as_used(aParentEntityIndex);
            mEntityRequestInfo(mRequestCounter, 0) = aParentEntityIndex;
            if (mChildEntityRank == EntityRank::NODE)
            {
                Integer* tIndex = mEntityTracker.get_index_pointer(aParentEntityIndex);
                Integer* tId = mEntityTracker.get_id_pointer(aParentEntityIndex);
                mPendingNodes(mRequestCounter).set_pending_node_info(tIndex, tId, aChildCoordsGlb, aParentTopology, aChildCoordsLoc);
            }
            mRequestCounter++;
        }

        // Get pointer to where the index will be put
        return mEntityTracker.get_index_pointer(aParentEntityIndex);
    }

    /*
     * Same as above but for when a secondary entity index is needed
     * Returns a pointer to the entity index
     */
    Integer*
    set_request_info(Integer aParentEntityIndex,
                     Integer aSecondaryEntityIdentifier,
                     Topology<Real,Integer,Real_Matrix,Integer_Matrix> const & aParentTopology,
                     moris::Matrix<Real,Real_Matrix>      const & aChildCoordsGlb,
                     moris::Matrix<Real,Real_Matrix>      const & aChildCoordsLoc,
                     moris::Matrix<Real,Real_Matrix>      const & aSensitivityDxDp = moris::Matrix<Real,Real_Matrix>(0,0),
                     moris::Matrix<Integer, Integer_Matrix> const & aNodeADVIndices  = moris::Matrix<Integer, Integer_Matrix>(0,0),
                     bool aHasDxdp = false,
                     bool aHasSparseDxDp = false)
    {
        // Ask entity tracker if this parent entity index has been used yet
        // For this function the flag is embedded in pIndId as Cell 0
        Cell<Integer*> pIdInd = mEntityTracker.is_parent_entity_used(aParentEntityIndex, aSecondaryEntityIdentifier);
        XTK_ASSERT(mNumChildrenAllowed != 1, "Coordinates submitted to set_request_info need to be one row (x,y,z)");

        if (pIdInd(0) == NULL)
        {
            if(mRequestCounter>=mEntityRequestInfo.n_rows())
            {
                XTK_ERROR<<"Warning: Not enough space allocated by constructor of Request handler, check aNumExpectedRequests. Dynamic Allocation"<<std::endl;
                Integer tNumRows  = mEntityRequestInfo.n_rows();
                Integer tNumCols  = mEntityRequestInfo.n_cols();
                mEntityRequestInfo.resize(2*tNumRows, tNumCols);
                mPendingNodes.resize(2*tNumRows, Pending_Node<Real, Integer,Real_Matrix, Integer_Matrix>());
            }


            mEntityRequestInfo(mRequestCounter, 0) = aParentEntityIndex;
            mEntityRequestInfo(mRequestCounter, 1) = aSecondaryEntityIdentifier;
            if (mChildEntityRank == EntityRank::NODE)
            {
                mPendingNodes(mRequestCounter).set_pending_node_info(pIdInd(2), pIdInd(1), aChildCoordsGlb, aParentTopology, aChildCoordsLoc);

                if(aHasDxdp)
                {
                    mPendingNodes(mRequestCounter).set_sensitivity_dx_dp(aSensitivityDxDp);
                }

                if(aHasSparseDxDp)
                {
                    mPendingNodes(mRequestCounter).set_node_adv_indices(aNodeADVIndices);
                }
            }
            mRequestCounter++;
        }

        // Get pointer to where the index will be put
        return pIdInd(2);
    }

    void handle_requests(bool const & aCoordFlag,
                         bool const & mSameMesh,
                         bool const & aInterfaceFlag,
                         Geometry_Engine<Real, Integer, Real_Matrix, Integer_Matrix> & aGeometryEngine)
    {
        // Resize out the unused space
        mEntityRequestInfo.resize(mRequestCounter, 2);
        mPendingNodes.resize(mRequestCounter, Pending_Node<Real, Integer,Real_Matrix, Integer_Matrix>());

        // Initialize Active Process Managers
        int tProcSize = 0;
        int tProcRank = 0;
        MPI_Comm_size(get_comm(), &tProcSize);
        MPI_Comm_rank(get_comm(), &tProcRank);
        Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> tActiveSendProcs(true,mNumChildrenAllowed,tProcSize,mParentEntityRank,mMesh);
        Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> tActiveRecvProcs(false,mNumChildrenAllowed,tProcSize,mParentEntityRank,mMesh);

        // Sort Requests and assign
        this->sort_entity_requests_and_assign_locally_controlled_entity_information(aCoordFlag,tActiveSendProcs,tActiveRecvProcs);

        // tell mesh to communicate requests (mesh internally places ids and index in the entity tracker
        this->communicate_entity_requests(aCoordFlag,tActiveSendProcs,tActiveRecvProcs);

        // Batch create the new entities in the background mesh(commits pending children entities to external entities in mesh)
        mMesh.batch_create_new_nodes(mPendingNodes);

        aGeometryEngine.associate_new_nodes_with_geometry_object(mPendingNodes,aInterfaceFlag);


    }

    void
    mark_pending_nodes_as_interface_nodes(XTK_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & aXTKMesh,
                                          Integer aGeometryIndex)
    {
        for(Integer i = 0; i <mPendingNodes.size(); i++)
        {
            aXTKMesh.mark_node_as_interface_node(mPendingNodes(i).get_node_index(),aGeometryIndex);

            // Determine if this is an interface node for any of the previous geometries
            Topology<Real,Integer,Real_Matrix,Integer_Matrix> const & tParentTopo = mPendingNodes(i).get_parent_topology();

            moris::Matrix<Integer, Integer_Matrix> const & tParentNodesInds = tParentTopo.get_node_indices();
            for(Integer iG = 0; iG<aGeometryIndex; iG++)
            {
                // If both nodes are created on an interface, then this node is an interface node with respect to the same geometry
                bool tIsInterfaceNode = true;
                for(Integer iN = 0; iN<tParentNodesInds.n_cols(); iN++)
                {
                    if(!aXTKMesh.is_interface_node(tParentNodesInds(0,iN),iG))
                    {
                        tIsInterfaceNode = false;
                        break;
                    }
                }

                if(tIsInterfaceNode)
                {
                    aXTKMesh.mark_node_as_interface_node(mPendingNodes(i).get_node_index(),iG);
                }
            }
        }


    }

    Integer get_num_requests() const
    {
        return mRequestCounter;
    }

    // -------------------DEBUG FUNCTIONS --------------------------
//    void
//    print();

private:
    Integer mRequestCounter;
    Integer mNumChildrenAllowed;
    Entity_Tracker<Real, Integer,Real_Matrix, Integer_Matrix> mEntityTracker;
    enum EntityRank mParentEntityRank;
    enum EntityRank mChildEntityRank;
    moris::Matrix<Integer, Integer_Matrix> mEntityRequestInfo;
    Cell<Pending_Node<Real, Integer,Real_Matrix, Integer_Matrix>> mPendingNodes;
    mesh::Mesh_Data<Real, Integer, Real_Matrix, Integer_Matrix> & mMesh;
    Cut_Mesh<Real, Integer, Real_Matrix, Integer_Matrix> & mXTKMesh;

    // Private Member function
private:

    /**
     * This function sorts
     */
    void sort_entity_requests_and_assign_locally_controlled_entity_information(bool aCoordFlag,
                                                                               Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> & aActiveSendProcs,
                                                                               Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> & aActiveRecvProcs)
    {
        int tProcSize = 0;
        int tProcRank = 0;

        MPI_Comm_size(get_comm(), &tProcSize);
        MPI_Comm_rank(get_comm(), &tProcRank);

        Integer tNumReqs = this->get_num_requests();

        // Allocate Glb entity Ids to each processor (1st MPI communication)
        // This currently assigns a block of continuous ids. This could be changed to return a list of non-contiguous Ids if id waste is excessive
        // the above change would require minimal changes to code on this side and slight changes in allocate_entity_ids
        Integer tLocalIdOffset = mMesh.allocate_entity_ids(tNumReqs, mChildEntityRank);
        Integer tLocalIndex = mMesh.get_first_available_index(mChildEntityRank);


        // Just do it in serial
        if (tProcSize == 1)
        {
            for (Integer i = 0; i < tNumReqs; i++)
            {

                // Assign a global Id and the local index
                mEntityTracker.set_child_entity_glb_id(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIdOffset);
                mEntityTracker.set_child_entity_lcl_index(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIndex);
                tLocalIdOffset++;
                tLocalIndex++;
                continue;
            }
        }

        // Go through the parallel process (which is more involved than the serial version because it requires knowledge of ownership)
        else
        {
            // Loop over requests and populate/ communicate assigned Ids
            for (Integer i = 0; i < tNumReqs; i++)
            {
                moris::Matrix<Integer, Integer_Matrix> tSharedProcs(1,1);
                mMesh.get_processors_whom_share_entity(mEntityRequestInfo(i, 0),mParentEntityRank,tSharedProcs);
//
                // Entity is not shared (these types of entities do not require any communication)
                if (((int)tSharedProcs(0, 0) == tProcRank) && (tSharedProcs.n_cols() == 1))
                {
                    // Assign a global Id and the local index
                    mEntityTracker.set_child_entity_glb_id(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIdOffset);
                    mEntityTracker.set_child_entity_lcl_index(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIndex);
                    tLocalIdOffset++;
                    tLocalIndex++;
                    continue;
                }
                // entity is shared (this requires a communication)
                else
                {
                    // Ask who owns the entity
                    int tOwnerProcRank = (int) mMesh.get_entity_parallel_owner_rank(mEntityRequestInfo(i, 0), mParentEntityRank);

                    // If the current processor owns the entity
                    if (tOwnerProcRank == tProcRank)
                    {
                        // Assign a global Id and index
                        mEntityTracker.set_child_entity_glb_id(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIdOffset);
                        mEntityTracker.set_child_entity_lcl_index(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIndex);

                        // Add entity Id to sending communication list
                        // It is important to send the ID because STK mesh can map the ID to a processor local index but a local index is basically garbage for other processors
                        for (Integer pr = 0; pr < tSharedProcs.n_cols(); pr++)
                        {
                            aActiveSendProcs.set_communication_info(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIdOffset, tSharedProcs(0, pr));
                        }

                        // advance local id and local index
                        tLocalIdOffset++;
                        tLocalIndex++;
                        continue;
                    }
                    // If the current processor does not own the entity
                                          

                    else if (tOwnerProcRank != tProcRank)
                        {
                            mEntityTracker.set_child_entity_glb_id(mEntityRequestInfo(i, 0),    mEntityRequestInfo(i, 1), tLocalIdOffset);
                            mEntityTracker.set_child_entity_lcl_index(mEntityRequestInfo(i, 0), mEntityRequestInfo(i, 1), tLocalIndex);
                            aActiveRecvProcs.set_communication_info(mEntityRequestInfo(i, 0),   mEntityRequestInfo(i, 1), tLocalIdOffset, tSharedProcs(0, 0));
                            tLocalIndex++;
                            tLocalIdOffset++;
                            continue;
                        }
                    
                }
            }
        }

        // Update the first available Inds
        aActiveSendProcs.condense_info();
        aActiveRecvProcs.condense_info();

        mMesh.update_first_available_index(tLocalIndex, mChildEntityRank);
    }
    //TODO: [MPI] This is only finished in serial
    void communicate_entity_requests(bool aCoordFlag,
                                     Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> & aActiveSendProcs,
                                     Active_Process_Manager<Real, Integer, Real_Matrix, Integer_Matrix> & aActiveRecvProcs)
    {

        //TODO: Come up with a designated MPI TAG method
        int tTag = 48;
        int tActiveProcessRank = 0;
        size_t tNumRows        = 0;
        size_t tNumColumns     = 0;
        Integer tParentIndex   = 0;
        Integer tSecondaryId   = 0;
        Integer tGlobalId      = 0;

        // Generate send and receive tags (could be changed to MPI_ANYTAG where the proc rank is the only decider)
        if(aActiveSendProcs.has_information())
        {
            Integer tNumSend = aActiveSendProcs.get_num_active_processors();

            for (Integer s = 0; s < tNumSend; s++)
            {

                // Get message information and size
                moris::Matrix<Integer, Integer_Matrix> & tSendMessage = aActiveSendProcs.get_comm_info(s);
                tNumRows = tSendMessage.n_rows();
                tNumColumns = tSendMessage.n_cols();

                // Get active Process to send message to
                tActiveProcessRank = aActiveSendProcs.get_active_processor_rank(s);

                // Send the message
                nonblocking_send(tSendMessage,tNumRows,tNumColumns,tActiveProcessRank,tTag);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);

        if(aActiveRecvProcs.has_information())
        {
            Integer tNumRecv = aActiveRecvProcs.get_num_active_processors();

            for (Integer r = 0; r < tNumRecv; r++)
            {
                // Get message information and size
                moris::Matrix<Integer, Integer_Matrix> & tRecvMessage = aActiveRecvProcs.get_comm_info(r);
                tNumRows = tRecvMessage.n_rows();

                // Get active Process to send message to
                tActiveProcessRank = aActiveRecvProcs.get_active_processor_rank(r);

                receive(tRecvMessage, tNumRows, tActiveProcessRank,tTag);

                tNumColumns = tRecvMessage.n_cols();
                for(Integer j = 0; j<tNumColumns; j++)
                {
                    tParentIndex = mMesh.get_loc_entity_index_from_entity_glb_id(tRecvMessage(0,j),mParentEntityRank);

                    tSecondaryId = tRecvMessage(1,j);
                    tGlobalId = tRecvMessage(2,j);

                    if(tParentIndex!=4294967295)
                    {
                    mEntityTracker.set_child_entity_glb_id(tParentIndex,tSecondaryId,tGlobalId);
                    }
                    // This is one type of hanging node. The processor did not expect to get a node here but did.
                }
            }
        }
        MPI_Barrier(get_comm());
    }

};
}




#endif /* SRC_XTK_CL_XTK_REQUEST_HANDLER_HPP_ */
