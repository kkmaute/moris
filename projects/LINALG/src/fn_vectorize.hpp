/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_vectorize.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_VECTORIZE_HPP_
#define PROJECTS_LINALG_SRC_FN_VECTORIZE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_vectorize_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_vectorize_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Generate a flattened version of matrix to either column (default) or row vector.
     *
     * for column vector: elements are copied column-wise, resulting in a column vector;
     *                    equivalent to concatenating all the columns of matrix
     *
     * for row vector:    elements are copied row-wise, resulting in a row vector;
     *                    equivalent to concatenating all the rows of matrix
     *
     * @param[in] aA            Matrix.
     * @param[in] aToRowVector  flag to create row vector
     *
     * @return  row or column vector
     */

    template< typename Matrix_Type >
    auto
    vectorize(
            const moris::Matrix< Matrix_Type > & aA,
            bool                                 aToRowVector = false)
    -> decltype( vectorize( aA.matrix_data(), aToRowVector) )
    {
        return vectorize( aA.matrix_data(), aToRowVector );
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_VECTORIZE_HPP_ */

