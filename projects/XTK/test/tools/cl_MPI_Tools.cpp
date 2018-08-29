/*
 * cl_MPI_Tools.cpp
 *
 *  Created on: Jun 27, 2017
 *      Author: ktdoble
 */

#include <mpi.h>
#include <vector>
#include "catch.hpp"

#include "tools/cl_MPI_Tools.hpp"

#include <containers/cl_XTK_Cell.hpp>
#include "core/xtk_typedefs.hpp"
#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "linalg/cl_XTK_Matrix.hpp"

#include "linalg_typedefs.hpp"




TEST_CASE("Gather","[MPI][GATHER][n2]")
{
    /*
     * If these tests fail its because there is probably an incorrect data type in their MPI call
     */
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    // Gather an Array
    if (tProcSize > 1)
    {
        SECTION("Gather of uint")
        {
            xtk::Cell<xtk::uint> tResultOfGather;
            xtk::Cell<xtk::uint> tMessage(1);
            int tMessageSize = 1;

            if (tProcRank == 0)
            {
                tMessage(0) = 15;
            }

            else if (tProcRank == 1)
            {
                tMessage(0) = 3;
            }

            else
            {
                tMessage(0) = 1;
            }

            xtk::gather(tMessage, tResultOfGather);

            if (tProcRank == 0)
            {
                REQUIRE(tResultOfGather(0) == 15);
                REQUIRE(tResultOfGather(1) == 3);
            }
        }

        SECTION("Gather of size_t")
        {
//            xtk::Cell<size_t> tResultOfGather;
//            xtk::Cell<size_t> tMessage;
//            int tMessageSize = 1;
//
//            if (tProcRank == 0)
//            {
//                tMessage(0) = 15;
//            }
//
//            else if (tProcRank == 1)
//            {
//                tMessage(0) = 3;
//            }
//
//            else
//            {
//                tMessage(0) = 1;
//            }
//
//            xtk::gather(tMessage, tResultOfGather);
//
//            if (tProcRank == 0)
//            {
//                REQUIRE(tResultOfGather(0) == 15);
//                REQUIRE(tResultOfGather(1) == 3);
//            }
        }
    }
}

TEST_CASE("Scatter","[MPI][SCATTER]")
{
    /*
     * If these tests fail its because there is probably an incorrect data type in their MPI call
     */
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    xtk::Cell<xtk::uint> tResultOfScatter(1);
    xtk::Cell<xtk::uint> tMessage(tProcSize);

    xtk::uint tVal = 24;
    xtk::uint tVal0 = tVal;
    if(tProcRank == 0)
    {
        for(xtk::uint i = 0; i<(xtk::uint)tProcSize; i++)
        {
            tMessage(i) = tVal;
            tVal++;
        }
    }

    xtk::scatter(tMessage, tResultOfScatter);
    REQUIRE(tResultOfScatter(0) == (xtk::uint)tProcRank+tVal0);
}


TEST_CASE("Send and Receive of XTK Matrix Class","[MPI][SEND][RECEIVE][n2]")
{
    int tProcRank = 0;
    int tProcSize = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &tProcRank);
    MPI_Comm_size(MPI_COMM_WORLD, &tProcSize);

    if(tProcSize==2)
    {
        // Initialize Matrix Manager
        moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tMatrix1(1,1);
        moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tMatrix2(1,1);
        if(tProcRank==0)
        {
            // Small Matrix Communication
            tMatrix1 = moris::Mat_New<size_t, xtk::Default_Matrix_Integer>(
                    {{10,11,12},
                {13,14,15},
                {22,26,36}});

            xtk::nonblocking_send(tMatrix1,3,3,1,0);

            // Large Matrix Communication Test
            xtk::nonblocking_send(tMatrix2,1000,1000,1,1);
        }

        else if (tProcRank == 1)
        {
            moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tMatrixRec1(1, 1, 0);
            moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tExpectedMatrixRec1(
                    {{ 10, 11, 12},
                { 13, 14, 15},
                { 22, 26, 36}});
            xtk::receive(tMatrixRec1, 3, 0, 0);
            CHECK(xtk::equal_to(tMatrixRec1,tExpectedMatrixRec1));

            moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tMatrixRec2(1, 1, 0);
            moris::Mat_New<size_t, xtk::Default_Matrix_Integer> tExpectedRecMatrix2(1000,1000,23);
            xtk::receive(tMatrixRec2,1000,0,1);
            CHECK(xtk::equal_to(tMatrixRec2,tExpectedRecMatrix2));
        }

        // So nothing goes out of scope before receiving
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

