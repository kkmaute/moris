/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_det.hpp
 *
 */

#pragma once

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_det_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_trace_Arma.hpp"
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Calculate the determinant of square matrix.
     *
     *@param[in] aA A given square matrix
     *
     * Example:
     * @include LNA/src/fn_det.inc
     *
     */
    template< typename Matrix_Type >
    auto
    trace( Matrix< Matrix_Type > const & aA )
    -> decltype( trace(aA.matrix_data()) ) const
    {
        return trace(aA.matrix_data());
    }

    template< typename Matrix_Type >
    auto
    trace( Matrix< Matrix_Type > & aA )
    -> decltype( trace(aA.matrix_data()) )
    {
        return trace(aA.matrix_data());
    }
}

