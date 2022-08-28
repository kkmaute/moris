/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * LINALG_Tutorial.cpp
 *
 */

#include "catch.hpp"

/*
 * The core moris typedefs
 */
#include "typedefs.hpp"
/*
 * The matrix include allows for use or moris::Matrix class
 */
#include "cl_Matrix.hpp"

/*
 * The linalg type include allows for the use of various matrix typedefs which
 * are defined based on a chose third party library (TPL). Note: the developer
 * does not need to be aware of the TPL being used.
 * i.e. A dense dynamic real matrix typedef is as follows:
 *      moris::DDRMat
 * where dense indicates there is no sparsity structure, dynamic means
 * there are no specified number of colums and rows are compile time,
 * real is the matrix data type.
 */
#include "linalg_typedefs.hpp"

/*
 * All linear algebra functions that are used in a file need to be included explicitly
 */
#include "fn_print.hpp"
#include "op_times.hpp"
#include "fn_eye.hpp"
#include "fn_trans.hpp"

/*
 * These includes are for using the std timing
 */
#include <iostream>
#include <ctime>

namespace moris
{
TEST_CASE("Basics of Linear Algebra Tutorial",
          "[LINALG_TUTORIAL]")
{
    // Since this test prints to console and does not explicitly test any
    // functions, it can be suppressed with this flag.
    bool tTutorialOn = false;

    if(tTutorialOn)
    {
        /*!
         * \section Introduction
         * The linear algebra package in MORIS is designed to enable the developer's of
         * MORIS to program similar to the syntax of MATLAB Linear algebra and to provide a
         * suite of operations/functions. Rather than write these operations functions,
         * we interface with third party libraries (TPLs) to provide the implementation.
         */

        /*!
         * \section Initialization
         *  A matrix is can be intitialized in a couple of ways.
         *  - 1.) with the number of columns and number of rows
         *  - 2.) with the number of columns and number of rows and fill value
         *  - 3.) with a nested std::initializer list
         */

        /*!
         * This allocates a 3x3 matrix but provides no values
         * \code{.cpp}
         * Matrix<DDRMat> tMat1(3,3);
         * \endcode
         */
        Matrix<DDRMat> tMat1(3,3);

        /*!
         * the print free function accepts a matrix and a title
         * and prints the contents of the matrix to the console
         */
        print(tMat1,"tMat1");

        /*!
         *  This allocates a 3x3 matrix and fills it with 1.0\n
         *  [ 1.0, 1.0, 1.0 \n
         *    1.0, 1.0, 1.0 \n
         *    1.0, 1.0, 1.0];
         * \code{.cpp}
         * Matrix<DDRMat> tMat2(3,3,0.0);
         * \endcode
         */
        Matrix<DDRMat> tMat2(3,3,0.0);
        print(tMat2,"tMat2");

        /*!
         * This allocates a 3x3 matrix \n
         *  [ 1.0, 2.0, 3.0 \n
         *    4.0, 5.0, 6.0 \n
         *    7.0, 8.0, 9.0];
         * \code{.cpp}
         * Matrix<DDRMat> tMat3({{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0}});
         * \endcode
         */
        Matrix<DDRMat> tMat3({{1.0,2.0,3.0},{4.0,5.0,6.0},{7.0,8.0,9.0}});
        print(tMat3,"tMat3");

        /*!
         * \section Free Functions and Operators
         * The linear algebra class features a suite free functions and operators.
         * Each free function and operator needs to be defined such that is handles all possible
         * combination between a moris Matrix and library expression templates.
         * Some examples of free functions are the identity function and the
         * transpose of a matrix.
         * \code{.cpp}
         * Matrix<DDRMat> tEye;
         * eye(3,3,tEye);
         * Matrix<DDRMat> tTransMat3 = trans(tMat3);
         * \endcode
         */
        Matrix<DDRMat> tEye;
        eye(3,3,tEye);
        print(tEye," A 3x3 identity matrix");

        Matrix<DDRMat> tTransMat3 = trans(tMat3);
        print(tTransMat3," Transpose of tMat3");

        /*!
         * The linalg package defines various operators between two matrices
         * i.e.
         * To use a multiply operator on two matrices
         * \code{.cpp}
         * Matrix<DDRMat> tMat4 = tMat2*tMat3;
         * \endcode
         */
        Matrix<DDRMat> tMat4 = tMat2*tMat3;
        print(tMat4, "tMat4");

        /*!
         * For a list of all operators and free functions, the reader is referred
         * to the LINALG directory where all headers starting with fn_ are free functions
         * and all headers starting with op_ are operators on matrices. Additionally,
         * test cases in LINALG/test show each of the operators/free functions in use
         * and their expected outputs.
         */

        /*!
         * \section On Expression Templating
         * Operations, like the times operator, and most free functions, like the
         * transpose operator do not explicitly return a matrix, both Eigen and Armadillo
         * TPLs leverage an expression template framework. The details of expression
         * templating will not be discussed here, but in general expression templates
         * allow for compiler optimization of series a linear algebra operations and prevent
         * unnecessary copies of matrixes to be made along the way. To show the benefit of
         * expression templates a comparison of times with and without using
         * expression templates is shown here. for the operations
         */

        size_t tNumIts = 100000;

        /*
         * No expression template loop
         */
        std::clock_t tStart;
        tStart = std::clock();
        Matrix<DDRMat> tMat6(3,3);
        for(size_t i = 0; i<tNumIts; i++)
        {
            Matrix<DDRMat> tMat5 = tEye*tMat3;
            Matrix<DDRMat> tMat2Trans = trans(tMat2);
            tMat6 = tMat5*tMat2Trans;
        }
        std::cout << "Time w/o expression templating: " << (std::clock() - tStart) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

        /*
         * With Expression templating
         */
        tStart = std::clock();
        for(size_t i = 0; i<tNumIts; i++)
        {
            tMat6 = tEye*tMat3*trans(tMat2);
        }
        std::cout << "Time w/ expression templating: " << (std::clock() - tStart) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl;

    }
}
}

