/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MPI_Tools.hpp
 *
 */

#ifndef SRC_TOOLS_CL_MPI_TOOLS_HPP_
#define SRC_TOOLS_CL_MPI_TOOLS_HPP_

#include <stdlib.h>
#include <mpi.h>

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"

#include "cl_Communication_Tools.hpp"

namespace moris
{
    /* ---------------------------------------------------------------------
     * Gather vectors of type T of all processors on root processor (default: 0)
     * Note: vectors need to have same size
     *
     * @param aBuffer input vector
     * @param aResult output vector (only valid on root processor)
     */
    template <typename T>
    inline
    void gather(
            Vector<T> & aBuffer,
            Vector<T> & aResult,
            int              aRootProc = 0)
    {
        // get processor rank and size
        int tProcRank = 0;
        int tProcSize = 0;

        MPI_Comm_rank(moris::get_comm(), &tProcRank);
        MPI_Comm_size(moris::get_comm(), &tProcSize);

        // get data type of vector
        MPI_Datatype tType = get_comm_datatype ( ( T ) 0 );

        // allocate proper size for root processor
        int tSizeofBuffer = aBuffer.size();
        if (tProcRank == aRootProc)
        {
            aResult.resize(tProcSize * tSizeofBuffer, 0);
        }

        // check for equal size of vectors across all processors
        MORIS_ASSERT( max_all( tSizeofBuffer ) == tSizeofBuffer,
                "gather - vector size differs across processors.\n");

        // mpi call
        MPI_Gather(
                aBuffer.data().data(),
                tSizeofBuffer,
                tType,
                aResult.data().data(),
                tSizeofBuffer,
                tType,
                aRootProc,
                moris::get_comm());
    }

    /* ---------------------------------------------------------------------
     * Scatter vectors of type T form root processor  (default: 0) to all processors
     * Note: vectors need to have same size
     *
     * @param aBuffer input vector (only valid on root processor)
     * @param aResult output vector (needs to be size correctly externally)
     */
    template <typename T>
    inline
    void scatter(
            Vector<T> & aBuffer,
            Vector<T> & aResult,
            int              aRootProc = 0)
    {
        // get processor rank and size
        int tProcRank = 0;
        int tProcSize = 0;

        MPI_Comm_rank(moris::get_comm(), &tProcRank);
        MPI_Comm_size(moris::get_comm(), &tProcSize);

        // get data type of vector
        MPI_Datatype tType = get_comm_datatype ( ( T ) 0 );

        // get size of receive buffer (needs to be set externally)
        int tSizeofBuffer = aResult.size();

        // check for equal size of cells across all processors
        MORIS_ASSERT( max_all( tSizeofBuffer ) == tSizeofBuffer,
                "scatter - cell size differs across processors.\n");

        // check for correct size of cell to be scattered
        MORIS_ASSERT( tProcRank == aRootProc ? aBuffer.size() == (uint)(tProcSize * tSizeofBuffer) : true,
                "scatter - inconsistent size of sending and receiving cells.");

        MPI_Scatter(
                aBuffer.data().data(),
                tSizeofBuffer,
                tType,
                aResult.data().data(),
                tSizeofBuffer,
                tType,
                aRootProc,
                moris::get_comm());
    }

    // ---------------------------------------------------------------------
    // FIXME: the following routine need proper commenting and error checking
    // ---------------------------------------------------------------------

    /**
     * Requirement: Column major matrix
     */
    template <typename Size_T_Matrix>
    MPI_Request nonblocking_send(
            moris::Matrix<Size_T_Matrix> const & aSendingMatrix,
            size_t                               aNumRows,
            size_t                               aNumColumns,
            int                                  aReceivingProc,
            int                                  aTag)
    {
        // get number of entries in matrix
        int tNumToSend = aSendingMatrix.numel();

        // initialize request
        MPI_Request tRequest;

        // non-blocking send
        MPI_Isend(
                aSendingMatrix.data(),
                tNumToSend,
                moris::get_comm_datatype(aSendingMatrix(0,0)),
                aReceivingProc,
                aTag,
                moris::get_comm(),
                &tRequest);

        return tRequest;
    }

    //------------------------------------------------------------------------------

    template <typename Size_T_Matrix>
    inline
    void receive(
            moris::Matrix<Size_T_Matrix> & aReceivingMatrix,
            size_t aNumRows,
            int aSendingProc,
            int aTag)
    {
        // Check that number of rows is not zero
        MORIS_ERROR( aNumRows > 0,
               "receive(moris::Matrix<T>) - number of rows cannot be zero." );

        // Check and wait until message is available
        MPI_Status tStatus;
        MPI_Probe(aSendingProc, aTag, moris::get_comm(), &tStatus);

        // Get size of message
        int tNumSent = 0;
        MPI_Get_count(
                &tStatus,
                moris::get_comm_datatype(aReceivingMatrix(0,0)),
                &tNumSent);

        // Check that size of matrix is larger zero
        MORIS_ERROR( tNumSent > 0,
                "receive(moris::Matrix<T>) - size of matrix to be received is zero." );

        // Compute number of columns
        size_t tNumColumns = tNumSent / aNumRows;

        // Check for proper number of columns
        MORIS_ERROR( (int)(tNumColumns*aNumRows) == tNumSent,
               "receive(moris::Matrix<T>) - inconsistent matrix sizes." );

        // Resize the matrix
        aReceivingMatrix.resize(aNumRows, tNumColumns);

        // Receive message
        MPI_Recv(
                aReceivingMatrix.data(),
                tNumSent,
                moris::get_comm_datatype(aReceivingMatrix(0,0)),
                aSendingProc,
                aTag,
                moris::get_comm(),
                &tStatus);
    }

    //------------------------------------------------------------------------------

}

#endif /* SRC_TOOLS_CL_MPI_TOOLS_HPP_ */

