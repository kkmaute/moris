#include <catch.hpp>

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
        std::string tFilePath = "matrixtest_" + std::to_string( moris::par_rank() ) + ".dat";

        // seed random number generator using time
        std::srand( time ( NULL ) );

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
