/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_join_vert.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_
#define PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_

// MORIS library header files.
#include <type_traits>
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_join_vert_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_join_vert_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief vertical concatenation
     *
     * @param[in] aA The First Matrix
     * @param[in] aB The Second Matrix
     *
     * @return join the corresponding columns of the given matrices
     *
     * @note The given matrices must have the same number of columns
     */
    template< typename Matrix_Type >
    auto
    join_vert(
            Matrix< Matrix_Type > const & aA,
            Matrix< Matrix_Type > const & aB )
            -> decltype( join_vert( aA.matrix_data(), aB.matrix_data() ) )
    {
        MORIS_ASSERT( aA.n_cols() == 0 || aB.n_cols() == 0 || aA.n_cols() == aB.n_cols(),
                "moris::join_vert - The number of columns in the two matrices must be equal." );

        return join_vert( aA.matrix_data(), aB.matrix_data() );
    }
    // A native concatenation expression may borrow its operands. Materialize
    // when either wrapper is temporary so returning the result remains safe.
    // The evaluated native type also grows fixed-size inputs to the joined size.
    template< typename Matrix_Type >
    auto join_vert( Matrix< Matrix_Type >&& aA, const Matrix< Matrix_Type >& aB )
    {
        auto tJoined      = join_vert( static_cast< const Matrix< Matrix_Type >& >( aA ), aB );
        using Result_Type = std::decay_t< decltype( tJoined.eval() ) >;
        return Matrix< Result_Type >( tJoined );
    }

    template< typename Matrix_Type >
    auto join_vert( const Matrix< Matrix_Type >& aA, Matrix< Matrix_Type >&& aB )
    {
        auto tJoined      = join_vert( aA, static_cast< const Matrix< Matrix_Type >& >( aB ) );
        using Result_Type = std::decay_t< decltype( tJoined.eval() ) >;
        return Matrix< Result_Type >( tJoined );
    }

    template< typename Matrix_Type >
    auto join_vert( Matrix< Matrix_Type >&& aA, Matrix< Matrix_Type >&& aB )
    {
        auto tJoined      = join_vert( static_cast< const Matrix< Matrix_Type >& >( aA ),
                static_cast< const Matrix< Matrix_Type >& >( aB ) );
        using Result_Type = std::decay_t< decltype( tJoined.eval() ) >;
        return Matrix< Result_Type >( tJoined );
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_JOIN_VERT_HPP_ */
