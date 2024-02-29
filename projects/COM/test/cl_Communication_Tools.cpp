/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Communication_Tools.cpp
 *
 */

#include <catch.hpp>

#include "moris_typedefs.hpp"    // COR/src
#include "banner.hpp"      // COR/src

#include "cl_Communication_Tools.hpp"      // COM/src
#include "cl_Communication_Manager.hpp"    // COM/src

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"
#include "fn_norm.hpp"

// Sam included headers for the sleep timer
#include <iostream>    // std::cout, std::endl

namespace moris
{
    TEST_CASE( "moris::send_mat_to_proc",
            "[comm]" )
    {

        /*SECTION( "moris::Mat send and receive test" )
    {
        // a vector to be communicated
        Matrix<DDUMat> tVec;

        // a matrix to be communicated
        Matrix< DDRMat > tMat;

        // task for first proc
        if ( moris::par_rank() == 0 )
        {
            // create vector to send
            tVec.resize( 6, 1 );
            tVec( 0 ) =  4;
            tVec( 1 ) =  8;
            tVec( 2 ) = 15;
            tVec( 3 ) = 16;
            tVec( 4 ) = 23;
            tVec( 5 ) = 42;

            // send vector to last proc
            send_mat_to_proc( tVec, moris::par_size() - 1 );

            // receive matrix from last proc
            recv_mat_from_proc( tMat, moris::par_size() - 1 );
        }

        /// task for last proc
        if ( moris::par_rank() == moris::par_size() - 1 )
        {

            std::cout<<" last"<<std::endl;
            // create matrix to send
            tMat.resize( 3, 2 );
            tMat( 0, 0 ) =  3;
            tMat( 1, 0 ) =  5;
            tMat( 2, 0 ) =  7;
            tMat( 0, 1 ) = 11;
            tMat( 1, 1 ) = 13;
            tMat( 2, 1 ) = 17;

            // create expected solution for vector to receive
            Matrix<DDUMat> tSolution(6, 1);
            tSolution( 0 ) =  4;
            tSolution( 1 ) =  8;
            tSolution( 2 ) = 15;
            tSolution( 3 ) = 16;
            tSolution( 4 ) = 23;
            tSolution( 5 ) = 42;

            send_mat_to_proc( tMat, 0 );
            recv_mat_from_proc( tVec, 0 );

            // calculate error vector
            Matrix<DDUMat> tError = tVec - tSolution;

            // make sure that received vector is correct
            REQUIRE( tError.max() == 0 && tError.min() == 0 );
        }

        // wait for all procs to be done
        moris::barrier();

        // check if matrix is correct
        if ( moris::par_rank() == 0 )
        {
            // create expected solution for matrix to receive
            Matrix< DDRMat > tSolution( 3, 2 );
            tSolution( 0, 0 ) =  3;
            tSolution( 1, 0 ) =  5;
            tSolution( 2, 0 ) =  7;
            tSolution( 0, 1 ) = 11;
            tSolution( 1, 1 ) = 13;
            tSolution( 2, 1 ) = 17;

            // calculate error matrix
            Matrix< DDRMat > tError = tMat - tSolution;

            // make sure that received matrix is correct
            REQUIRE( tError.max() == 0 && tError.min() == 0 );
        }

    } */
    }

    TEST_CASE( "moris::sum_min_max",
            "[comm],[allreduce]" )
    {
        // check min, max, sum of real value
        real tLocalReal = (real)par_rank();

        real tGlobalRealMax = max_all( tLocalReal );
        real tGlobalRealMin = min_all( tLocalReal );
        real tGlobalRealSum = sum_all( tLocalReal );

        real tCheckRealSum = 0.0;
        for ( uint i = 0; i < (uint)par_size(); ++i ) tCheckRealSum += (real)i;

        REQUIRE( tGlobalRealMax == (real)( par_size() - 1 ) );
        REQUIRE( tGlobalRealMin == 0.0 );
        REQUIRE( tGlobalRealSum == tCheckRealSum );

        // check sum of real matrix
        Matrix< DDRMat > tLocalRealMat( 10, 5, (real)par_rank() );
        Matrix< DDRMat > tCheckRealMat( 10, 5, tCheckRealSum );

        Matrix< DDRMat > tGlobalRealMatSum = sum_all_matrix( tLocalRealMat );

        REQUIRE( norm( tGlobalRealMatSum - tCheckRealMat ) < 1e-12 );
    }

    TEST_CASE( "moris::proc_cart",
            "[comm],[proc_cart]" )
    {

        if ( moris::par_size() == 1 or moris::par_size() == 2 or moris::par_size() == 4 )
        {
            uint             tNumberOfDimensions = 3;
            Matrix< DDUMat > tProcDims;
            Matrix< DDUMat > tProcDimsCompare;
            Matrix< DDUMat > tError;
            Matrix< DDUMat > tProcCoords;
            Matrix< IdMat >  tProcNeighbors;

            for ( uint iDecompMethod = 0; iDecompMethod <= 2; ++iDecompMethod )
            {
                if ( moris::par_size() == 1 )
                {
                    tProcDims        = { { 1 }, { 1 }, { 1 } };
                    tProcDimsCompare = tProcDims;
                }
                else if ( moris::par_size() == 2 )
                {
                    if ( iDecompMethod == 1 )
                    {
                        tProcDimsCompare = { { 2 }, { 1 }, { 1 } };
                    }
                    else if ( iDecompMethod == 2 )
                    {
                        tProcDims        = { { 1 }, { 2 }, { 1 } };
                        tProcDimsCompare = tProcDims;
                    }
                    else if ( iDecompMethod == 0 )
                    {
                        tProcDims        = { { 1 }, { 1 }, { 2 } };
                        tProcDimsCompare = tProcDims;
                    }
                }
                else
                {
                    if ( iDecompMethod == 1 )
                    {
                        tProcDimsCompare = { { 2 }, { 2 }, { 1 } };
                    }
                    else if ( iDecompMethod == 2 )
                    {
                        tProcDims        = { { 1 }, { 2 }, { 2 } };
                        tProcDimsCompare = tProcDims;
                    }
                    else if ( iDecompMethod == 0 )
                    {
                        tProcDims        = { { 2 }, { 1 }, { 2 } };
                        tProcDimsCompare = tProcDims;
                    }
                }

                create_proc_cart(
                        iDecompMethod,
                        tNumberOfDimensions,
                        tProcDims,
                        tProcCoords,
                        tProcNeighbors );

                tError = tProcDims - tProcDimsCompare;
                REQUIRE( ( tError.max() == 0 && tError.min() == 0 ) );
            }
        }
    }
}    // namespace moris
