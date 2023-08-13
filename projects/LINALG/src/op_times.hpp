/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * op_times.hpp
 *
 */

#ifndef PROJECTS_LINALG_OP_TIMES_HPP_
#define PROJECTS_LINALG_OP_TIMES_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/op_times_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/op_times_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Matrix multiplication
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function returns a matrix such that @f$ \mathbf{C}_{ij} = \mathbf{A}_{ik} \mathbf{B}_{kj} @f$
     * where C is the product of A and B. If A is an m-by-p and B is a p-by-n matrix, then C is an m-by-n matrix.
     *
     * Example:
     * @include LNA/src/op_times.inc
     *
     */
    template< typename Matrix_Type_A, typename Matrix_Type_B >
    auto
    operator*( Matrix< Matrix_Type_A >& aA,
            Matrix< Matrix_Type_B >&    aB )
            -> decltype( aA.matrix_data() * aB.matrix_data() )
    {

        MORIS_ASSERT( aA.n_cols() == aB.n_rows(), "Dimension mismatch in matrix multiplication. %-5zu vs %-5zu", aA.n_cols(), aB.n_rows() );
        return aA.matrix_data() * aB.matrix_data();
    }

    template< typename Matrix_Type >
    auto
    operator*( typename Matrix< Matrix_Type >::Data_Type& aA,
            Matrix< Matrix_Type >&                        aB )
            -> decltype( scalar_times( aA, aB.matrix_data() ) )
    {
        return scalar_times( aA, aB.matrix_data() );
    }

    template< typename Matrix_Type_A, typename Matrix_Type_B >
    auto
    operator*( const Matrix< Matrix_Type_A >& aA,
            const Matrix< Matrix_Type_B >&    aB )
            -> decltype( aA.matrix_data() * aB.matrix_data() )
    {
        MORIS_ASSERT( aA.n_cols() == aB.n_rows(), "Dimension mismatch in matrix multiplication. %-5zu vs %-5zu", aA.n_cols(), aB.n_rows() );
        return aA.matrix_data() * aB.matrix_data();
    }

    template< typename Matrix_Type >
    auto
    operator*( const typename Matrix< Matrix_Type >::Data_Type& aA,
            const Matrix< Matrix_Type >&                        aB )
            -> decltype( scalar_times( aA, aB.matrix_data() ) )
    {
        return scalar_times( aA, aB.matrix_data() );
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_OP_TIMES_HPP_ */
