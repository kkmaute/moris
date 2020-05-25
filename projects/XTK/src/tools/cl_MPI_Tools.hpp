/*
 * cl_MPI_Tools.hpp
 *
 *  Created on: Jun 27, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOOLS_CL_MPI_TOOLS_HPP_
#define SRC_TOOLS_CL_MPI_TOOLS_HPP_

#include <stdlib.h>
#include <mpi.h>

#include "xtk_typedefs.hpp"
#include "cl_Cell.hpp"

#include "cl_Matrix.hpp"

#include "cl_Communication_Tools.hpp"


namespace xtk
{

inline
void gather(moris::Cell<xtk::uint> & aBuffer, moris::Cell<xtk::uint> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(moris::get_comm(), &tProcRank);
    MPI_Comm_size(moris::get_comm(), &tProcSize);

    int tSizeofBuffer = aBuffer.size();
    if (tProcRank == 0)
    {
        aResult.resize(tProcSize * tSizeofBuffer, 0);
    }

    MPI_Gather(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED, tRoot, moris::get_comm());
}

inline
void gather(moris::Cell<moris::moris_id> & aBuffer,
            moris::Cell<moris::moris_id> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(moris::get_comm(), &tProcRank);
    MPI_Comm_size(moris::get_comm(), &tProcSize);

    int tSizeofBuffer = aBuffer.size();
    if (tProcRank == 0)
    {
        aResult.resize(tProcSize * tSizeofBuffer, 0);
    }

    MPI_Datatype tDataType = moris::get_comm_datatype(aBuffer(0));

    MPI_Gather(aBuffer.data().data(), tSizeofBuffer, tDataType, aResult.data().data(), tSizeofBuffer, tDataType, tRoot, moris::get_comm());
}

/**
 * Require the result to be sized appropriately
 */
inline
void scatter(moris::Cell<xtk::uint> & aBuffer, moris::Cell<xtk::uint> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(moris::get_comm(), &tProcRank);
    MPI_Comm_size(moris::get_comm(), &tProcSize);

//    MORIS_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED, tRoot, moris::get_comm());
}

inline
void scatter(moris::Cell<moris::moris_id> & aBuffer,
             moris::Cell<moris::moris_id> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(moris::get_comm(), &tProcRank);
    MPI_Comm_size(moris::get_comm(), &tProcSize);

//    MORIS_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Datatype tDataType = moris::get_comm_datatype(aBuffer(0));

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, tDataType, aResult.data().data(), tSizeofBuffer, tDataType, tRoot, moris::get_comm());
}



/**
 * Requirement: Column major matrix
 */
template <typename Size_T_Matrix>
void nonblocking_send(moris::Matrix<Size_T_Matrix> const & aSendingMatrix,
                      size_t aNumRows,
                      size_t aNumColumns,
                      int aReceivingProc,
                      int aTag)
{
    int tNumToSend = aSendingMatrix.numel();

    MPI_Request tRequest;


    std::cout<<"NONBLOCKING SEND FROM "<<par_rank()<<"TO "<<aReceivingProc<<" WITH TAG "<<aTag<<std::endl;

    MPI_Isend(aSendingMatrix.data(), tNumToSend, moris::get_comm_datatype(aSendingMatrix(0,0)), aReceivingProc, aTag, moris::get_comm(), &tRequest);
}

inline
bool
sent_message_exists(int aOtherProc,
                    int aTag)
{

    int flag = 1000;
    MPI_Status tStatus;
    MPI_Iprobe(aOtherProc, aTag, moris::get_comm(), &flag, &tStatus);

    std::cout<<"Flag = "<<flag<< " | aOtherProc = "<<aOtherProc<<"| my rank = "<<moris::par_rank()<<" | Tag = "<<aTag<<std::endl;

    if(flag)
    {
        return true;
    }
    else
    {
        return false;
    }
}

template <typename Size_T_Matrix>
inline
void receive(moris::Matrix<Size_T_Matrix> & aReceivingMatrix,
             size_t aNumRows,
             int aSendingProc,
             int aTag)
{
    MORIS_ERROR(sent_message_exists(aSendingProc,aTag),"Trying to receive a message that does not exists");
    MPI_Status tStatus;
    int tNumSent = 0;
    MPI_Get_count(&tStatus, moris::get_comm_datatype(aReceivingMatrix(0,0)), &tNumSent);


    size_t tNumColumns = tNumSent / aNumRows;

    // Resize the matrix
    aReceivingMatrix.resize(aNumRows, tNumColumns);

    MPI_Recv(aReceivingMatrix.data(), tNumSent, moris::get_comm_datatype(aReceivingMatrix(0,0)), aSendingProc, aTag, moris::get_comm(), &tStatus);

}

template <typename Size_T_Matrix>
inline
void receive_col_known(moris::Matrix<Size_T_Matrix> & aReceivingMatrix,
             size_t aNumCols,
             int aSendingProc,
             int aTag)
{

    MORIS_ERROR(sent_message_exists(aSendingProc,aTag),"Trying to receive a message that does not exists");
    MPI_Status tStatus;
    int tNumSent = 0;
    MPI_Get_count(&tStatus, moris::get_comm_datatype(aReceivingMatrix(0,0)), &tNumSent);

    size_t tNumRows = tNumSent / aNumCols;

    // Resize the matrix
    aReceivingMatrix.resize(tNumRows, aNumCols);

    MPI_Recv(aReceivingMatrix.data(), tNumSent, moris::get_comm_datatype(aReceivingMatrix(0,0)), aSendingProc, aTag, moris::get_comm(), &tStatus);

}

inline
moris::moris_id
allocate_ids(moris::size_t aNumIdstoAllocate,
             moris::size_t tStart = 1)
{
    int tProcRank = moris::par_rank();
    int tProcSize = moris::par_size();

    // size_t is defined as uint here because of aNumRequested
    //Initialize gathered information outputs (information which will be scattered across processors)
    moris::Cell<moris::moris_id> aGatheredInfo;
    moris::Cell<moris::moris_id> tFirstId(1);
    moris::Cell<moris::moris_id> tNumIdsRequested(1);

    tNumIdsRequested(0) = (moris::moris_id)aNumIdstoAllocate;

    xtk::gather(tNumIdsRequested,aGatheredInfo);

    moris::Cell<moris::moris_id> tProcFirstID(tProcSize);

    moris::moris_id tFirstAvailableId = tStart;

    if(tProcRank == 0)
    {
        // Loop over entities print the number of entities requested by each processor
        for (int iProc = 0; iProc < tProcSize; ++iProc)
        {
            // Give each processor their desired amount of IDs
            tProcFirstID(iProc) = tFirstAvailableId;

            // Increment the first available node ID
            tFirstAvailableId = tFirstAvailableId+aGatheredInfo(iProc);
        }

    }

    xtk::scatter(tProcFirstID,tFirstId);


    return tFirstId(0);
}




}


#endif /* SRC_TOOLS_CL_MPI_TOOLS_HPP_ */
