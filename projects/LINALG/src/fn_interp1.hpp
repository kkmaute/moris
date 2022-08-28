/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_interp1.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_INTERP1_HPP_
#define PROJECTS_LINALG_SRC_FN_INTERP1_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

// function does not exist in Eigen

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_interp1_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief 1D data interpolation. Given a 1D function specified in vectors X and Y (where X specifies locations and Y
     * specifies the corresponding values), generate vector YI which contains interpolated values at locations XI
     *
     * @param[in] aX A vector of 1D locations
     * @param[in] aY A vector of Values at the X locations
     * @param[in] aXI A vector of Values of 1D locations to interpolate at
     * @param[out] aYI A vector of interpolated values
     *
     */
    template< typename Matrix_Type >
    auto
    interp1( Matrix< Matrix_Type > const & aX,
             Matrix< Matrix_Type > const & aY,
             Matrix< Matrix_Type > const & aXI,
             Matrix< Matrix_Type >  & aYI )
    -> decltype( interp1(aX.matrix_data(),aY.matrix_data(),aXI.matrix_data(),aYI.matrix_data())  )
    {
        return interp1(aX.matrix_data(),aY.matrix_data(),aXI.matrix_data(),aYI.matrix_data());
    }
}
#endif /* PROJECTS_LINALG_SRC_FN_INTERP1_HPP_ */

