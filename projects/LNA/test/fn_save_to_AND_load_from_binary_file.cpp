#include <string>


// Third-party header files.
#include <catch.hpp>

// MORIS project header files.
#include "typedefs.hpp" // COR/src
#include "cl_Communication_Tools.hpp" // COM/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_load_vector_from_binary_file.hpp" // LNA/src
#include "fn_save_vector_to_binary_file.hpp" // LNA/src

TEST_CASE(
        "moris::save_to_binary_file",
        "[linalgebra],[save_to_binary_file],[load_from_binary_file]" )
{
    SECTION( "testing binary IO of vector" )
    {

        uint tNumberOfSamples = 750;

        // determine file path
        std::string tFilePath = "matrixtest_" + std::to_string( moris::par_rank() ) + ".dat";

        // seed random number generator using time
        std::srand( time ( NULL ) );

        // assign vector
        moris::Mat< moris:: real > tOutVec( tNumberOfSamples, 1 );

        // write random values into vector
        for( moris::uint k=0; k<tNumberOfSamples; ++k )
        {
            tOutVec( k ) = std::rand();
        }

        // save vector into file
        moris::save_vector_to_binary_file( tFilePath , tOutVec ) ;

        // input vector
        moris::Mat< moris:: real > tInVec;

        // load vector from file
        moris::load_vector_from_binary_file ( tFilePath , tInVec ) ;

        // assert that both vectors have the same length
        REQUIRE( tInVec.length() == tOutVec.length());

        // calculate error
        moris::real tError = 0;
        for ( uint k=0; k<tNumberOfSamples; ++k )
        {
            tError += std::pow( tOutVec( k ) - tInVec( k ), 2 );
        }

        // create norm
        tError = std::pow( tError, 0.5 );

        // make sure that norm is less than an epsilon
        REQUIRE( tError <= 1E-6 );

    }
}
