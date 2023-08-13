/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isvector.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_ISVECTOR_HPP_
#define PROJECTS_LINALG_SRC_FN_ISVECTOR_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_isvector_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_isvector_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Checks if a given object is vector
 *
 * @param[in] aA Input object to be checked
 *
 * This function operates on given object (i.e. matrix, vector) and returns
 * true if the object can be interpreted as a vector (either column or row
 * vector). It also return false if the object isn't vector.

 *
 */
    template< typename Matrix_Type >
    auto
    isvector( Matrix< Matrix_Type > const & aA )
    -> decltype( isvector( aA.matrix_data() ) )
    {
        return isvector( aA.matrix_data() );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_ISVECTOR_HPP_ */

