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
#include "banner.hpp"            // COR/src

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

    //------------------------------------------------------------------------------

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

    //------------------------------------------------------------------------------

    TEST_CASE( "moris::proc_cart",
            "[comm],[proc_cart]" )
    {

        if ( moris::par_size() == 1 or moris::par_size() == 2 or moris::par_size() == 4 )
        {
            uint             tNumberOfDimensions = 3;
            Vector< uint >   tProcDims;
            Vector< uint > tProcDimsCompare;
            Matrix< DDUMat > tError;
            Matrix< DDUMat > tProcCoords;
            Matrix< IdMat >  tProcNeighbors;

            for ( uint iDecompMethod = 0; iDecompMethod <= 2; ++iDecompMethod )
            {
                if ( moris::par_size() == 1 )
                {
                    tProcDims        = { 1, 1, 1 };
                    tProcDimsCompare = tProcDims;
                }
                else if ( moris::par_size() == 2 )
                {
                    if ( iDecompMethod == 1 )
                    {
                        tProcDimsCompare = { 2, 1, 1 };
                    }
                    else if ( iDecompMethod == 2 )
                    {
                        tProcDims        = { 1, 2, 1 };
                        tProcDimsCompare = tProcDims;
                    }
                    else if ( iDecompMethod == 0 )
                    {
                        tProcDims        = { 1, 1, 2 };
                        tProcDimsCompare = tProcDims;
                    }
                }
                else
                {
                    if ( iDecompMethod == 1 )
                    {
                        tProcDimsCompare = { 2, 2, 1 };
                    }
                    else if ( iDecompMethod == 2 )
                    {
                        tProcDims        = { 1, 2, 2 };
                        tProcDimsCompare = tProcDims;
                    }
                    else if ( iDecompMethod == 0 )
                    {
                        tProcDims        = { 2, 1, 2 };
                        tProcDimsCompare = tProcDims;
                    }
                }

                create_proc_cart(
                        iDecompMethod,
                        tNumberOfDimensions,
                        tProcDims,
                        tProcCoords,
                        tProcNeighbors );

                for ( uint iCheckIndex = 0; iCheckIndex < tProcDims.size(); iCheckIndex++ )
                {
                    REQUIRE( tProcDims( iCheckIndex ) == tProcDimsCompare( iCheckIndex ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "moris::gatherv_mats",
            "[comm],[gather_mats]" )
    {
        // a matrix to be communicated
        Matrix< DDRMat > tMat;

        // vector of matrices to be gathered
        Vector< Matrix< DDRMat > > tAllMats;

        // task for first proc
        if ( moris::par_rank() == 0 )
        {
            // create vector to send
            tMat = { { 11, 12 }, { 21, 22 }, { 31, 32 } };
        }
        else
        {
            tMat.set_size( 2, par_rank(), (real)par_rank() );
        }

        gatherv_mats( tMat, tAllMats );

        // make sure that received matrix is correct
        if ( moris::par_rank() == 0 )
        {
            // check for size of vector of matrices
            REQUIRE( tAllMats.size() == (uint)par_size() );

            // check first matrix
            REQUIRE( norm( tMat - tAllMats( 0 ) ) < 1e-12 );

            // check last matrix
            for ( sint i = 1; i < par_size(); ++i )
            {
                Matrix< DDRMat > tOffProcMat( 2, i, (real)i );

                REQUIRE( norm( tOffProcMat - tAllMats( i ) ) < 1e-12 );
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "moris::scatter_mats",
            "[comm],[scatter_mats]" )
    {
        // a matrix to be communicated
        Matrix< DDRMat > tMat;

        Vector< Matrix< DDRMat > > tAllMats;

        if ( moris::par_rank() == 0 )
        {
            // vector of matrices to be gathered
            tAllMats.resize( par_size() );

            // create vector to send
            tAllMats( 0 ) = { { 11, 12 }, { 21, 22 }, { 31, 32 } };

            for ( sint i = 1; i < par_size(); ++i )
            {
                tAllMats( i ).set_size( 2, i, (real)i );
            }
        }

        scatterv_mats( tAllMats, tMat );

        if ( moris::par_rank() == 0 )
        {
            Matrix< DDRMat > tReference = { { 11, 12 }, { 21, 22 }, { 31, 32 } };

            // check first matrix
            REQUIRE( norm( tMat - tReference ) < 1e-12 );
        }
        else
        {
            Matrix< DDRMat > tReference( 2, par_rank(), (real)par_rank() );

            REQUIRE( norm( tMat - tReference ) < 1e-12 );
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "moris::braodcast_mat",
            "[comm],[braodcast_mat]" )
    {
        // a matrix to be communicated
        Matrix< DDRMat > tMat;

        Matrix< DDRMat > tReference = { { 11, 12 }, { 21, 22 }, { 31, 32 } };

        if ( moris::par_rank() == 0 )
        {
            // create vector to send
            tMat = tReference;
        }

        broadcast_mat( tMat );

        REQUIRE( norm( tMat - tReference ) < 1e-12 );
    }

}    // namespace moris
