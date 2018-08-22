/*
 * cl_XTK_Matrix.cpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */

#include "linalg/cl_XTK_Matrix_Base.hpp"

#include "core/xtk_typedefs.hpp"

#ifdef MORIS_USE_EIGEN
#include "linalg/cl_XTK_Matrix_Eigen.hpp"
#else
#include "linalg/cl_XTK_Matrix_Teuchos.hpp"
#endif

#include "linalg/cl_XTK_Matrix_Base_Utilities.hpp"
#include "catch.hpp"

namespace xtk
{
TEST_CASE("XTK Linear Algebra Matrix Base Tests","[LINALG]")
{
    SECTION("Matrix Tests using Default definition"){
    // Create matrices
#ifdef MORIS_USE_EIGEN
    xtk::Matrix_Base_Eigen<real, Default_Matrix_Real> tMatrix1(1,2);
#else
    xtk::Matrix_Base_Teuchos<real, Default_Matrix_Real> tMatrix1(1,2);
#endif


    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix2 = tMatrix1.create();
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix3 = tMatrix1.create(5,4);
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix4 = tMatrix1.create(5,4,-1);

    // Check number of columns
    REQUIRE(tMatrix1.get_num_columns() == 2);
    REQUIRE(tMatrix2->get_num_columns() == 2);
    REQUIRE(tMatrix3->get_num_columns() == 4);

    // Check number of rows
    REQUIRE(tMatrix1.get_num_rows() == 1);
    REQUIRE(tMatrix2->get_num_rows() == 1);
    REQUIRE(tMatrix3->get_num_rows() == 5);

    // Check number of elements in matrices
    REQUIRE(tMatrix1.get_num_matrix_elements() == 2);
    REQUIRE(tMatrix2->get_num_matrix_elements() == 2);
    REQUIRE(tMatrix3->get_num_matrix_elements() == 20);

    // Resize tMatrix2 and check again
    tMatrix2->resize(12,4);
    REQUIRE(tMatrix2->get_num_rows() == 12);
    REQUIRE(tMatrix2->get_num_columns() == 4);
    REQUIRE(tMatrix2->get_num_matrix_elements() == 48);

    // Add Matrix Data by row index and location
    tMatrix1(0,0) = 1;
    tMatrix1(0,1) = 2;
    REQUIRE(tMatrix1(0,0) == 1);
    REQUIRE(tMatrix1(0,1) == 2);

    // Add Matrix data to std::shared_ptr<Matrix<type>>
    (*tMatrix2)(0,0)= 1;
    (*tMatrix2)(0,1)= 2;
    REQUIRE((*tMatrix2)(0,0) == 1);
    REQUIRE((*tMatrix2)(0,1) == 2);

    // Check filling operation of tMat4 (internal to create());
    REQUIRE((*tMatrix4)(0,0) == -1);

    // Check explicit fill call
    tMatrix4->fill(48);
    REQUIRE((*tMatrix4)(0,0) == 48);

    // Set and get Columns
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix5 = tMatrix1.create(4,4,10);
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrixRow1 = tMatrix1.create(1,4,0);
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrixColumn1 = tMatrix1.create(4,1,0);

    tMatrix5->get_row(2,(*tMatrixRow1));
    REQUIRE((*tMatrixRow1)(0,0) == 10);

    tMatrix5->get_column(1,*tMatrixColumn1);
    REQUIRE((*tMatrixColumn1)(1,0) == 10);

    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrixRow2 = tMatrix1.create(1,4,5);
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrixColumn2 = tMatrix1.create(4,1,7);

    tMatrix5->set_row(3,(*tMatrixRow2));
    tMatrix5->set_column(3,(*tMatrixColumn2));

    REQUIRE((*tMatrix5)(3,2) == 5);
    REQUIRE((*tMatrix5)(3,3) == 7);

    // Test maximum and minimum values
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix6 = tMatrix1.create(10,7,0);
    (*tMatrix6)(3,5) = 10;
    (*tMatrix6)(5,5) = -11;

    REQUIRE(tMatrix6->get_max_value() == 10);
    REQUIRE(tMatrix6->get_min_value() == -11);

    // Create a matrix using a standard initializer list

    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix7 = tMatrix1.create({{1,2,3},{4,5,6},{7,8,9}});

    REQUIRE((*tMatrix7)(0,0) = 1);
    REQUIRE((*tMatrix7)(1,1) = 5);
    REQUIRE((*tMatrix7)(2,2) = 9);

    // Check Data function
    const real* tMatrix7Data = tMatrix7->data();

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
    std::shared_ptr<Matrix_Base<real, Default_Matrix_Real>> tMatrix7Copy = tMatrix7->copy();

    CHECK(xtk::equal_to(*tMatrix7,*tMatrix7Copy));

    // Modify the original and see if the copy changes (it shouldn't)
    (*tMatrix7)(0,0) = 44;

    CHECK_FALSE(xtk::equal_to(*tMatrix7,*tMatrix7Copy));

    // Index out of bounds
    REQUIRE_THROWS((*tMatrix7)(3,3) = 0);

    // Initializer list with a mistake in it
    REQUIRE_THROWS(tMatrix1.create({{1,2,3,4},{4,5,6},{7,8,9}}));
}
}

}
