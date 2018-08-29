/*
 * cl_Matrix.cpp
 *
 *  Created on: Aug 27, 2018
 *      Author: doble
 */
#include <catch.hpp>
#include "fn_equal_to.hpp" //ALG

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris
{
TEST_CASE("MORIS Linear Algebra Matrix Tests","[MATRIX]")
{
    SECTION("Matrix Tests using default"){
        // Create matrix base
        Matrix< real, DDRMat > tMatrix1(1,2);
        Matrix< real, DDRMat > tMatrix2(0,0);
        Matrix< real, DDRMat > tMatrix3(5,4);
        Matrix< real, DDRMat > tMatrix4(5,4,-1);

        // Check number of columns
        REQUIRE(tMatrix1.n_cols() == 2);
        REQUIRE(tMatrix2.n_cols() == 0);
        REQUIRE(tMatrix3.n_cols() == 4);

        // Check number of rows
        REQUIRE(tMatrix1.n_rows() == 1);
        REQUIRE(tMatrix2.n_rows() == 0);
        REQUIRE(tMatrix3.n_rows() == 5);

        // Check number of elements in matrices
        REQUIRE(tMatrix1.numel() == 2);
        REQUIRE(tMatrix2.numel() == 0);
        REQUIRE(tMatrix3.numel() == 20);

        // Resize tMatrix2 and check again
        tMatrix2.resize(12,4);
        REQUIRE(tMatrix2.n_rows() == 12);
        REQUIRE(tMatrix2.n_cols() == 4);
        REQUIRE(tMatrix2.numel() == 48);

        // Add Matrix Data by row index and location
        tMatrix1(0,0) = 1;
        tMatrix1(0,1) = 2;
        REQUIRE(tMatrix1(0,0) == 1);
        REQUIRE(tMatrix1(0,1) == 2);

        // Add Matrix data to std::shared_ptr<Matrix<type>>
        tMatrix2(0,0)= 1;
        tMatrix2(0,1)= 2;
        REQUIRE(tMatrix2(0,0) == 1);
        REQUIRE(tMatrix2(0,1) == 2);

        // Check filling operation of tMat4 (internal to create());
        REQUIRE(tMatrix4(0,0) == -1);

        // Check explicit fill call
        tMatrix4.fill(48);
        REQUIRE(tMatrix4(0,0) == 48);

        // Set and get Columns
        Matrix< real, DDRMat >  tMatrix5(4,4,10);

        Matrix< real, DDRMat > tMatrixRow1 = tMatrix5.get_row(2);
        REQUIRE(tMatrixRow1(0,0) == 10);

        Matrix< real, DDRMat >  tMatrixColumn1 = tMatrix5.get_column(1);
        REQUIRE(tMatrixColumn1(1,0) == 10);

        Matrix< real, DDRMat > tMatrixRow2(1,4,5);
        Matrix< real, DDRMat > tMatrixColumn2(4,1,7);

        tMatrix5.set_row(3,tMatrixRow2);
        tMatrix5.set_column(3,tMatrixColumn2);

        REQUIRE(tMatrix5(3,2) == 5);
        REQUIRE(tMatrix5(3,3) == 7);

        // Test maximum and minimum values
        Matrix< real, DDRMat > tMatrix6(10,7,0);
        tMatrix6(3,5) = 10;
        tMatrix6(5,5) = -11;

        REQUIRE(tMatrix6.max() == 10);
        REQUIRE(tMatrix6.min() == -11);

        // Create a matrix using a standard initializer list

        Matrix< real, DDRMat > tMatrix7({{1,2,3},{4,5,6},{7,8,9}});

        REQUIRE(tMatrix7(0,0) = 1);
        REQUIRE(tMatrix7(1,1) = 5);
        REQUIRE(tMatrix7(2,2) = 9);

        // Check Data function
        const real* tMatrix7Data = tMatrix7.data();

        // Column Major Data Structure Check
        REQUIRE(tMatrix7Data[0] == 1);
        REQUIRE(tMatrix7Data[1] == 4);
        REQUIRE(tMatrix7Data[2] == 7);
        REQUIRE(tMatrix7Data[3] == 2);
        REQUIRE(tMatrix7Data[4] == 5);
        REQUIRE(tMatrix7Data[5] == 8);
        REQUIRE(tMatrix7Data[6] == 3);
        REQUIRE(tMatrix7Data[7] == 6);
        REQUIRE(tMatrix7Data[8] == 9);

        // Copying an existing matrix
        Matrix< real, DDRMat > tMatrix7Copy = tMatrix7.copy();

//        CHECK(xtk::equal_to(tMatrix7, tMatrix7Copy));

        // Modify the original and see if the copy changes (it shouldn't)
        tMatrix7(0,0) = 44;

//        CHECK_FALSE(xtk::equal_to(tMatrix7,tMatrix7Copy));
//
//        // Index out of bounds
//        REQUIRE_THROWS(tMatrix7(3,3) = 0);
//
//        // Initializer list with a mistake in it
//        REQUIRE_THROWS(Matrix< real, DDRMat >({{1,2,3,4},{4,5,6},{7,8,9}}));


    }
}
}

