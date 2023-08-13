/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_rank.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_RANK_HPP_
#define PROJECTS_LINALG_SRC_FN_RANK_HPP_

#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_rank_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_rank_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Return the rank of matrix
     *
     * @return Return the rank of matrix aA
     *
     * @note The computation is based on singular value decomposition;
     * if the decomposition fails, a std::runtime_error exception is thrown
     */
    template< typename Matrix_Type >
    auto
    rank(Matrix< Matrix_Type > const & aA)
    -> decltype( rank( aA.matrix_data()) )
    {
        return rank( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_RANK_HPP_ */

