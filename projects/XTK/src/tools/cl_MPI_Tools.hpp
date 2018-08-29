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

#include "core/xtk_typedefs.hpp"
#include "containers/cl_XTK_Cell.hpp"

#include "cl_Matrix.hpp"
#include "assert/fn_xtk_assert.hpp"



namespace xtk
{

MPI_Comm
inline get_comm()
{
  MPI_Comm tComm = MPI_COMM_WORLD;
  return tComm;
}
void set_comm(MPI_Comm aComm);

int
get_rank(MPI_Comm aComm)
{
    int tProcRank = 0;
    MPI_Comm_rank(aComm, &tProcRank);
    return tProcRank;
}

int
get_size(MPI_Comm aComm)
{
    int tProcSize = 0;
    MPI_Comm_size(aComm, &tProcSize);
    return tProcSize;
}

void gather(xtk::Cell<xtk::uint> & aBuffer, xtk::Cell<xtk::uint> & aResult)
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

void gather(xtk::Cell<xtk::size_t> & aBuffer, xtk::Cell<xtk::size_t> & aResult)
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

    MPI_Gather(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED_LONG_LONG, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED_LONG_LONG, tRoot, MPI_COMM_WORLD);
}

/**
 * Require the result to be sized appropriately
 */
void scatter(xtk::Cell<xtk::uint> & aBuffer, xtk::Cell<xtk::uint> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

//    XTK_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED, tRoot, MPI_COMM_WORLD);
}

void scatter(xtk::Cell<xtk::size_t> & aBuffer, xtk::Cell<xtk::size_t> & aResult)
{
    int tProcRank = 0;
    int tProcSize = 0;
    int tRoot = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

//    XTK_ASSERT(aResult.size() != 0, "aResult Buffer needs to externally allocate size");

    int tSizeofBuffer = aResult.size();

    MPI_Scatter(aBuffer.data().data(), tSizeofBuffer, MPI_UNSIGNED_LONG_LONG, aResult.data().data(), tSizeofBuffer, MPI_UNSIGNED_LONG_LONG, tRoot, MPI_COMM_WORLD);
}



/**
 * Requirement: Column major matrix
 */
template <typename Size_T_Matrix>
void nonblocking_send(moris::Mat_New<xtk::size_t, Size_T_Matrix> const & aSendingMatrix,
                      size_t aNumRows,
                      size_t aNumColumns,
                      int aReceivingProc,
                      int aTag)
{
    int tNumToSend = aNumRows * aNumColumns;

    MPI_Request tRequest;

    MPI_Isend(aSendingMatrix.data(), tNumToSend, MPI_UNSIGNED_LONG_LONG, aReceivingProc, aTag, get_comm(), &tRequest);
}

template <typename Size_T_Matrix>
void receive(moris::Mat_New<xtk::size_t, Size_T_Matrix> & aReceivingMatrix,
             size_t aNumRows,
             int aSendingProc,
             int aTag)
{
    // The number of rows is known but the number of columns depends on the sending message
    // Need to probe
    MPI_Status tStatus;
    int tFlag = 0;
    int tNumSent = 0;
    MPI_Iprobe(aSendingProc, aTag, get_comm(), &tFlag, &tStatus);

    if(tStatus.MPI_ERROR==0)
    {
        MPI_Get_count(&tStatus, MPI_UNSIGNED_LONG_LONG, &tNumSent);

        // Allocate Buffer space
        xtk::size_t* tBuffer = new xtk::size_t[tNumSent];

        size_t tNumColumns = tNumSent / aNumRows;

        // Resize the matrix
        aReceivingMatrix.resize(aNumRows, tNumColumns);

        MPI_Recv(tBuffer, tNumSent, MPI_UNSIGNED_LONG_LONG, aSendingProc, aTag, get_comm(), &tStatus);

        // Modify the provided matrix
        size_t k = 0;
        for (size_t c = 0; c < tNumColumns; c++)
        {
            for (size_t r = 0; r < aNumRows; r++)
            {
                aReceivingMatrix(r, c) = tBuffer[k];
                k++;
            }
        }

        delete [] tBuffer;

    }

    else
    {
        std::cout<<"No message received! Expected from processor "<<aSendingProc<< " to " << get_rank(get_comm())<<std::endl;
        aReceivingMatrix.resize(0,0);
    }
}
}


#endif /* SRC_TOOLS_CL_MPI_TOOLS_HPP_ */
