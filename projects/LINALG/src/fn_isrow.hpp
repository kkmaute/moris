/*
 * fn_isrow.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_FN_ISROW_HPP_
#define PROJECTS_LINALG_SRC_FN_ISROW_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_isrow_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_isrow_Arma.hpp"
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
