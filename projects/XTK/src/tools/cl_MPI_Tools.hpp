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
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    int tSizeofBuffer = aBuffer.size();
    if (tProcRank == 0)
    {
        aResult.resize(tProcSize * tSizeofBuffer, 0);
    }

    MPI_Gather(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED, tRoot, MPI_COMM_WORLD);
}

inline
void gather(moris::Cell<moris::moris_id> & aBuffer,
            moris::Cell<moris::moris_id> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    int tSizeofBuffer = aBuffer.size();
    if (tProcRank == 0)
    {
        aResult.resize(tProcSize * tSizeofBuffer, 0);
    }

    MPI_Datatype tDataType = moris::get_comm_datatype(aBuffer(0));

    MPI_Gather(aBuffer.data().data(), tSizeofBuffer, tDataType, aResult.data().data(), tSizeofBuffer, tDataType, tRoot, MPI_COMM_WORLD);
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
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

//    MORIS_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED, tRoot, MPI_COMM_WORLD);
}

inline
void scatter(moris::Cell<moris::moris_id> & aBuffer,
             moris::Cell<moris::moris_id> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

//    MORIS_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Datatype tDataType = moris::get_comm_datatype(aBuffer(0));

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, tDataType, aResult.data().data(), tSizeofBuffer, tDataType, tRoot, MPI_COMM_WORLD);
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


    MPI_Isend(aSendingMatrix.data(), tNumToSend, moris::get_comm_datatype(aSendingMatrix(0,0)), aReceivingProc, aTag, MPI_COMM_WORLD, &tRequest);
}

inline
bool
sent_message_exists(int aOtherProc,
                    int aTag)
{

    int flag = 1000;
    MPI_Status tStatus;
    MPI_Iprobe(aOtherProc, aTag, moris::get_comm(), &flag, &tStatus);


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
    // The number of rows is known but the number of columns depends on the sending message
    // Need to probe
    MPI_Status tStatus;
    int tFlag = 0;
    int tNumSent = 0;
    MPI_Iprobe(aSendingProc, aTag, MPI_COMM_WORLD, &tFlag, &tStatus);

    if(tStatus.MPI_ERROR==0)
    {
        MPI_Get_count(&tStatus, moris::get_comm_datatype(aReceivingMatrix(0,0)), &tNumSent);


        size_t tNumColumns = tNumSent / aNumRows;

        // Resize the matrix
        aReceivingMatrix.resize(aNumRows, tNumColumns);

        MPI_Recv(aReceivingMatrix.data(), tNumSent, moris::get_comm_datatype(aReceivingMatrix(0,0)), aSendingProc, aTag, MPI_COMM_WORLD, &tStatus);

    }

    else
    {
        std::cout<<"No message received! Expected from processor "<<aSendingProc<< " to " << moris::par_rank()<<std::endl;
        aReceivingMatrix.resize(0,0);
    }
}



}


#endif /* SRC_TOOLS_CL_MPI_TOOLS_HPP_ */
