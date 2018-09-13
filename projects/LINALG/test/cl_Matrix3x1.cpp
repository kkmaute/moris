/*
 * cl_Matrix3x1.cpp
 *
 *  Created on: Sep 11, 2018
 *      Author: schmidt
 */
#include <catch.hpp>
#include "fn_equal_to.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
TEST_CASE("MORIS 3x1 Matrix Test","[linalgebra],[MATRIX3x1]")
{
    SECTION("Matrix Tests using default"){
        // Create matrix base
        Matrix< real, F31RMat > tMatrix1( 3,1,1.1 );
        Matrix< uint, F31UMat > tMatrix2( 3,1,1 );

        // Check number of columns
        REQUIRE( tMatrix1.n_cols() == 1 );
        REQUIRE( tMatrix2.n_cols() == 1 );

        // Check number of rows
        REQUIRE( tMatrix1.n_rows() == 3 );
        REQUIRE( tMatrix2.n_rows() == 3 );

        // Check number of elements in matrices
        REQUIRE( tMatrix1.numel() == 3 );
        REQUIRE( tMatrix2.numel() == 3 );

        // Add Matrix Data by row index and location
        tMatrix1( 0, 0 ) = 3;
        tMatrix1( 1, 0 ) = 2;
        REQUIRE( tMatrix1( 0, 0 ) == 3 );
        REQUIRE( tMatrix1( 1, 0 ) == 2 );
        REQUIRE( tMatrix1( 2, 0 ) == 1.1 );

        // Add Matrix data to std::shared_ptr<Matrix<type>>
        tMatrix2( 0, 0 ) = 5;
        tMatrix2( 1, 0 ) = 2;
        REQUIRE( tMatrix2( 0, 0 ) == 5 );
        REQUIRE( tMatrix2( 1, 0 ) == 2 );
        REQUIRE( tMatrix2( 2, 0 ) == 1 );

        // Check explicit fill call
        tMatrix2.fill( 48 );
        REQUIRE( tMatrix2( 0, 0 ) == 48 );

        // Test maximum and minimum values
        Matrix< real, F31RMat > tMatrix3( 3,1,0 );
        tMatrix3( 1, 0 ) = 10;
        tMatrix3( 2, 0 ) = -11;

        REQUIRE( tMatrix3.max() == 10);
        REQUIRE( tMatrix3.min() == -11);

        // Create a matrix using a standard initializer list
        Matrix< real, F31RMat > tMatrix4( {{1},{2},{3}} );
        REQUIRE( tMatrix4( 0,0 ) = 1 );
        REQUIRE( tMatrix4( 1,0 ) = 2 );
        REQUIRE( tMatrix4( 2,0 ) = 3 );

        // Check Data function
        const real* tMatrix4Data = tMatrix4.data();

        // Column Major Data Structure Check
        REQUIRE( tMatrix4Data[0] == 1 );
        REQUIRE( tMatrix4Data[1] == 2 );
        REQUIRE( tMatrix4Data[2] == 3 );


        // Copying an existing matrix
        Matrix< real, F31RMat > tMatrix4Copy = tMatrix4.copy();
        CHECK( equal_to( tMatrix4( 0, 0 ), tMatrix4Copy( 0, 0 ) ) );

        // Modify the original and see if the copy changes (it shouldn't)
        tMatrix4( 0, 0 ) = 44;

        CHECK_FALSE( equal_to( tMatrix4( 0, 0 ), tMatrix4Copy( 0, 0 ) ) );

        // Index out of bounds
        REQUIRE_THROWS( tMatrix4( 3, 3 ) = 0 );

        // Initializer list with a mistake in it
//        REQUIRE_THROWS( (Matrix< real, F31RMat >( {{1,2},{4,5},{7,8}} )) );
//        REQUIRE_THROWS( (Matrix< real, F31RMat >( {{1},{4},{7},{7}} )) );


    }
}
}

