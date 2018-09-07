/*
 * fn_issquare.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */
#ifndef PROJECTS_LINALG_SRC_FN_ISSQUARE_HPP_
#define PROJECTS_LINALG_SRC_FN_ISSQUARE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_issquare_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_issquare_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Checks if a given object is square
 *
 * @param[in] aA Input object to be checked
 *
 * This function operates on given object (i.e. matrix, vector) and returns
 * true if the object can be interpreted as a square matrix (i.e. number of
 * columns equal to number of rows). It also return false if the object
 * isn't square matrix.
 *
 */
    template< typename Type, typename Matrix_Type >
    auto
    issquare( Matrix< Type, Matrix_Type  > const & aA )
    -> decltype( issquare( aA.matrix_data() ) )
    {
        return issquare( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_ISSQUARE_HPP_ */
