/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_eye.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_EYE_HPP_
#define PROJECTS_LINALG_SRC_FN_EYE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eye_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_eye_Arma.hpp"
#endif

namespace moris
{

    /**
     * @brief Writes Identity Matrix on given Matrix
     *
     * @param[in]  aNumRows Number of Rows.
     * @param[in]  aNumCols Number of Columns.
     * @param[out] aEyeMat  Identity matrix @f$ \mathbf{I}_{n}@f$
     * such that @f$ \mathbf{({I}_{n})}_{ij} =  \mathbf{\delta}_{ij} @f$\n
     */

    template< typename Matrix_Type >
    inline
    void
    eye(
            size_t const          & aNumRows,
            size_t const          & aNumCols,
            Matrix< Matrix_Type > & aEyeMat)
    {
        eye( aNumRows, aNumCols, aEyeMat.matrix_data() );
    }

    /**
     * @brief Creates an Identity Matrix
     *
     * @param[in] aNumRows Number of Rows.
     * @param[in] aNumCols Number of Columns.

     * @return Identity matrix @f$ \mathbf{I}_{n}@f$
     * such that @f$ \mathbf{({I}_{n})}_{ij} =  \mathbf{\delta}_{ij} @f$\n
     */

     inline
     auto
     eye(
             size_t const          & aNumRows,
             size_t const          & aNumCols)
     ->decltype ( linalg_internal::eye( aNumRows, aNumCols ) )
     {
         return linalg_internal::eye( aNumRows, aNumCols );
     }
}

#endif /* PROJECTS_LINALG_SRC_FN_EYE_HPP_ */

