/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_iscol.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_ISCOL_HPP_
#define PROJECTS_LINALG_SRC_FN_ISCOL_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_iscol_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_iscol_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Checks if matrix is a colvector, this means it has only one column.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Matrix_Type >
    auto
    iscol( Matrix< Matrix_Type > const & aA )
    -> decltype( iscol( aA.matrix_data() ) )
    {
        return iscol( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_ISCOL_HPP_ */

