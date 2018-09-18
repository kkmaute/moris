/*
 * fn_isempty.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_FN_ISEMPTY_HPP_
#define PROJECTS_LINALG_SRC_FN_ISEMPTY_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_isempty_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_isempty_Arma.hpp"
#endif


namespace moris
{
    /**
     * @brief Checks if a matrix is empty.
     *
     *@param[in] aA A given matrix
     *
     *
     */
    template< typename Matrix_Type >
    auto
    isempty( Matrix< Matrix_Type > const & aA )
    -> decltype( isempty( aA.matrix_data() ) )
    {
        return isempty( aA.matrix_data() );
    }
}



#endif /* PROJECTS_LINALG_SRC_FN_ISEMPTY_HPP_ */
