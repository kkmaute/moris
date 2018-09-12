/*
 * op_minus.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_OP_MINUS_HPP_
#define PROJECTS_LINALG_SRC_OP_MINUS_HPP_

#include "cl_Matrix.hpp"


#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_minus_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_minus_Arma.hpp"
#endif


namespace moris
{
    /**
     * @brief Subtraction of two arrays
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ij} - \mathbf{B}_{ij} @f$
     * where C is suntraction of B from A. If A is an m-by-p and B is a m-by-p matrix, then C is an m-by-p matrix.
     *
     * Example:
     * @include LNA/src/op_minus.inc
     *
     */
    template< typename Type, typename Matrix_Type >
    auto
    operator-(
            Matrix< Type, Matrix_Type > const & aA,
            Matrix< Type, Matrix_Type > const & aB )
    ->decltype( aA.matrix_data() - aB.matrix_data() )
    {
        return aA.matrix_data() - aB.matrix_data();
    }

    template< typename Type, typename Matrix_Type >
    auto
    operator-(
            Type const & aScalar,
            Matrix< Type, Matrix_Type > const & aMatrix )
    ->decltype( aScalar - aMatrix.matrix_data() )
    {
        return aScalar - aMatrix.matrix_data();
    }

    template< typename Type, typename Matrix_Type >
    auto
    operator-(
            Matrix< Type, Matrix_Type > const & aMatrix,
            Type const & aScalar
    )
    ->decltype( aMatrix.matrix_data() - aScalar )
    {
        return aMatrix.matrix_data() - aScalar;
    }
}

#endif /* PROJECTS_LINALG_SRC_OP_MINUS_HPP_ */
