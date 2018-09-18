/*
 * fn_eye.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_FN_EYE_HPP_
#define PROJECTS_LINALG_SRC_FN_EYE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_eye_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_eye_Arma.hpp"
#endif

namespace moris
{

/**
  * @brief Creates an Identity Matrix
  *
  * @param[in] aNumRows Number of Rows.
  * @param[in] aNumCols Number of Columns.
  * An identity matrix is generated when aNumRows = aNumCols
  *
  * @return Creates and identity matrix @f$ \mathbf{I}_{n}@f$
  * such that @f$ \mathbf{({I}_{n})}_{ij} =  \mathbf{\delta}_{ij} @f$\n
  * Generates a matrix with the elements along the main diagonal
  * set to one and off-diagonal elements set to zero.
  */

template< typename Matrix_Type >
inline
void
eye( size_t const &              aNumRows,
     size_t const &              aNumCols,
     Matrix< Matrix_Type > & aEyeMat)
 {
    eye( aNumRows, aNumCols, aEyeMat.matrix_data() );
 }

}


#endif /* PROJECTS_LINALG_SRC_FN_EYE_HPP_ */
