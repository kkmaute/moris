#include <catch.hpp>

#include "typedefs.hpp" // COR/src
#include "banner.hpp" // COR/src

#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Communication_Manager.hpp" // COM/src

#include "cl_Matrix.hpp" // LNA/src
#include "linalg_typedefs.hpp"
#include "op_minus.hpp"

namespace moris
{
TEST_CASE( "moris::send_mat_to_proc",
                "[comm]" )
{
    SECTION( "moris::Mat send and receive test" )
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

    }
}
}
