/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isrow.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_ISROW_HPP_
#define PROJECTS_LINALG_SRC_FN_ISROW_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_isrow_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_isrow_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Checks if matrix is a row vector, this means it has only one row.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Matrix_Type >
    auto
    isrow( Matrix< Matrix_Type > const & aA )
    -> decltype( isrow( aA.matrix_data() ) )
    {
        return isrow( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_ISROW_HPP_ */

