/*
 * fn_sort.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_FN_SORT_HPP_
#define PROJECTS_LINALG_SRC_FN_SORT_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_sort_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_sort_Arma.hpp"
#endif


namespace moris
{
/**
 * @brief Sorts colums of a vector/matrix.
 *
 * @param[in] aMat Vector.
 *
 * @return  sorted vector
 */
    template< typename Type, typename Matrix_Type >
    void
    sort( Matrix< Type, Matrix_Type  > const & aA,
          Matrix< Type, Matrix_Type  > & aSorted)
    {
        sort( aA.matrix_data(), aSorted );
    }

}

#endif /* PROJECTS_LINALG_SRC_FN_SORT_HPP_ */
