/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_save_and_load_matrix_binary_file.cpp
 *
 */

#include <catch.hpp>
#include <cstdio>

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Communication_Manager.hpp"
#include "cl_Communication_Tools.hpp"

#include "fn_save_matrix_to_binary_file.hpp"
#include "fn_load_matrix_from_binary_file.hpp"

extern moris::Comm_Manager gMorisComm;

TEST_CASE("moris::save_and_load_binary_file",
          "[linalgebra],[save_and_load_binary_file]")
{
        // set size of test matrix
        uint tNumberOfRows = 640;
        uint tNumberOfCols = 480;

        // determine file path
        std::string tFilePath = "matrixtest_" + std::to_string( moris::par_size() ) + "_" + std::to_string( moris::par_rank() ) + ".dat";

        // seed random number generator using time
        std::srand( time ( nullptr ) );

        // assign vector
        moris::Matrix< moris::DDRMat >  tOutMat( tNumberOfRows, tNumberOfCols );

        // write random values into vector
        for( uint j=0; j<tNumberOfCols; ++j )
        {
            for( uint i=0; i<tNumberOfRows; ++i )
            {
                tOutMat( i, j ) = std::rand();
            }
        }

        // save vector into file
        moris::save_matrix_to_binary_file( tOutMat, tFilePath ) ;

        // input vector
        moris::Matrix< moris::DDRMat > tInMat;

        // load vector from file
        moris::load_matrix_from_binary_file ( tInMat, tFilePath ) ;

        // assert that both vectors have the same rows and cols
        REQUIRE( tInMat.n_cols() == tOutMat.n_cols() );
        REQUIRE( tInMat.n_rows() == tOutMat.n_rows() );

        // calculate error
        moris::real tError = 0;

        for( uint j=0; j<tNumberOfCols; ++j )
        {
            for( uint i=0; i<tNumberOfRows; ++i )
            {
                tError += std::pow( tInMat( i, j ) - tOutMat( i, j ), 2 );
            }
        }

        // create norm
        tError = std::pow( tError, 0.5 );

        // make sure that norm is less than an epsilon
        REQUIRE( tError <= 1E-6 );
}


TEST_CASE( "LINALG binary IO crosses the full buffer boundary", "[linalgebra],[save_and_load_binary_file]" )
{
    // The implementation flushes 4 MiB at a time. Cover an exact flush and a
    // flush followed by a partial tail; the previous pointer bug affected both.
    for ( uint tTail : { 0u, 3u } )
    {
        const uint                     tCount = 4u * 1024u * 1024u / sizeof( moris::real ) + tTail;
        const std::string              tPath  = "matrix_buffer_" + std::to_string( moris::par_size() ) + "_" + std::to_string( moris::par_rank() ) + ".dat";
        moris::Matrix< moris::DDRMat > tWritten( tCount, 1 );
        for ( uint i = 0; i < tCount; ++i )
        {
            tWritten( i ) = static_cast< moris::real >( i ) / 8.0;
        }
        moris::save_matrix_to_binary_file( tWritten, tPath );
        moris::Matrix< moris::DDRMat > tRead;
        moris::load_matrix_from_binary_file( tRead, tPath );
        REQUIRE( tRead.n_rows() == tCount );
        REQUIRE( tRead.n_cols() == 1 );
        bool tEqual = true;
        for ( uint i = 0; i < tCount; ++i )
        {
            tEqual = tEqual && tRead( i ) == tWritten( i );
        }
        CHECK( tEqual );
        CHECK( std::remove( tPath.c_str() ) == 0 );
    }
}

TEST_CASE( "LINALG binary IO handles empty and complex matrices and open failures", "[linalgebra],[save_and_load_binary_file]" )
{
    const std::string                    tPath = "matrix_io_edges_" + std::to_string( moris::par_size() ) + "_" + std::to_string( moris::par_rank() ) + ".dat";
    const moris::Matrix< moris::DDRMat > tEmpty;
    moris::save_matrix_to_binary_file( tEmpty, tPath );
    moris::Matrix< moris::DDRMat > tRead( 2, 2, 1.0 );
    moris::load_matrix_from_binary_file( tRead, tPath );
    CHECK( tRead.n_rows() == 0 );
    CHECK( tRead.n_cols() == 0 );

    const moris::Matrix< moris::DDCMat > tComplex = { { moris::cplx( 1.0, 2.0 ), moris::cplx( 3.0, -4.0 ) } };
    moris::save_matrix_to_binary_file( tComplex, tPath );
    moris::Matrix< moris::DDCMat > tComplexRead;
    moris::load_matrix_from_binary_file( tComplexRead, tPath );
    REQUIRE( tComplexRead.n_rows() == 1 );
    REQUIRE( tComplexRead.n_cols() == 2 );
    CHECK( tComplexRead( 0 ) == moris::cplx( 1.0, 2.0 ) );
    CHECK( tComplexRead( 1 ) == moris::cplx( 3.0, -4.0 ) );
    REQUIRE( std::remove( tPath.c_str() ) == 0 );
    CHECK_THROWS_AS( moris::load_matrix_from_binary_file( tRead, tPath ), std::runtime_error );
    CHECK_THROWS_AS( moris::save_matrix_to_binary_file( tRead, tPath + "/missing.dat" ), std::runtime_error );
}
