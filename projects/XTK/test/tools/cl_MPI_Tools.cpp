/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MPI_Tools.cpp
 *
 */

#include <mpi.h>
#include <vector>
#include "catch.hpp"

#include <cl_XTK_Cell.hpp>
#include "xtk_typedefs.hpp"
#include "cl_XTK_Matrix_Base_Utilities.hpp"
#include "cl_Matrix.hpp"

#include "linalg_typedefs.hpp"

namespace moris::xtk
{

    TEST_CASE( "Gather", "[MPI][GATHER][n2]" )
    {
        /*
         * If these tests fail its because there is probably an incorrect data type in their MPI call
         */
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
        MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

        // Gather an Array
        if ( tProcSize > 1 )
        {
            SECTION( "Gather of uint" )
            {
                Vector< xtk::uint > tResultOfGather;
                Vector< xtk::uint > tMessage( 1 );
                int                 tMessageSize = 1;

                if ( tProcRank == 0 )
                {
                    tMessage( 0 ) = 15;
                }

                else if ( tProcRank == 1 )
                {
                    tMessage( 0 ) = 3;
                }

                else
                {
                    tMessage( 0 ) = 1;
                }

                xtk::gather( tMessage, tResultOfGather );

                if ( tProcRank == 0 )
                {
                    REQUIRE( tResultOfGather( 0 ) == 15 );
                    REQUIRE( tResultOfGather( 1 ) == 3 );
                }
            }

            SECTION( "Gather of size_t" )
            {
                //            Vector<size_t> tResultOfGather;
                //            Vector<size_t> tMessage;
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

    TEST_CASE( "Scatter", "[MPI][SCATTER]" )
    {
        /*
         * If these tests fail its because there is probably an incorrect data type in their MPI call
         */
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
        MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

        Vector< xtk::uint > tResultOfScatter( 1 );
        Vector< xtk::uint > tMessage( tProcSize );

        xtk::uint tVal  = 24;
        xtk::uint tVal0 = tVal;
        if ( tProcRank == 0 )
        {
            for ( xtk::uint i = 0; i < (xtk::uint)tProcSize; i++ )
            {
                tMessage( i ) = tVal;
                tVal++;
            }
        }

        xtk::scatter( tMessage, tResultOfScatter );
        REQUIRE( tResultOfScatter( 0 ) == (xtk::uint)tProcRank + tVal0 );
    }

    TEST_CASE( "Send and Receive of XTK Matrix Class", "[MPI][SEND][RECEIVE][n2]" )
    {
        int tProcRank = 0;
        int tProcSize = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &tProcRank );
        MPI_Comm_size( MPI_COMM_WORLD, &tProcSize );

        if ( tProcSize == 2 )
        {
            // Initialize Matrix Manager
            moris::Matrix< moris::DDSTMat > tMatrix1( 1, 1 );
            moris::Matrix< moris::DDSTMat > tMatrix2( 1, 1 );
            if ( tProcRank == 0 )
            {
                // Small Matrix Communication
                tMatrix1 = moris::Matrix< moris::DDSTMat >(
                        { { 10, 11, 12 },
                                { 13, 14, 15 },
                                { 22, 26, 36 } } );

                xtk::nonblocking_send( tMatrix1, 3, 3, 1, 0 );

                // Large Matrix Communication Test
                xtk::nonblocking_send( tMatrix2, 1000, 1000, 1, 1 );
            }

            else if ( tProcRank == 1 )
            {
                moris::Matrix< moris::DDSTMat > tMatrixRec1( 1, 1, 0 );
                moris::Matrix< moris::DDSTMat > tExpectedMatrixRec1(
                        { { 10, 11, 12 },
                                { 13, 14, 15 },
                                { 22, 26, 36 } } );
                xtk::receive( tMatrixRec1, 3, 0, 0 );
                CHECK( xtk::equal_to( tMatrixRec1, tExpectedMatrixRec1 ) );

                moris::Matrix< moris::DDSTMat > tMatrixRec2( 1, 1, 0 );
                moris::Matrix< moris::DDSTMat > tExpectedRecMatrix2( 1000, 1000, 23 );
                xtk::receive( tMatrixRec2, 1000, 0, 1 );
                CHECK( xtk::equal_to( tMatrixRec2, tExpectedRecMatrix2 ) );
            }

            // So nothing goes out of scope before receiving
            MPI_Barrier( MPI_COMM_WORLD );
        }
    }

}    // namespace moris::xtk
