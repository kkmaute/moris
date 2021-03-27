/*
 * fn_reshape.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_FN_RESHAPE_HPP_
#define PROJECTS_LINALG_SRC_FN_RESHAPE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_reshape_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_reshape_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief reshape vector.
     *
     * @param[in] aA     matrix
     * @param[in] aB     number of rows
     * @param[in] aC     number of columns
     *
     * @return  reshaped vector
     */

    template< typename Matrix_Type >
    auto
    reshape(
            const moris::Matrix< Matrix_Type > & aA,
            const size_t                       & aB,
            const size_t                       & aC)
    -> decltype( reshape( aA.matrix_data(),aB,aC) )
    {
        return reshape( aA.matrix_data(),aB,aC );
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_RESHAPE_HPP_ */
