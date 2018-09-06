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
#include "Arma_Impl/fn_reshape_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief find vector.
 *
 * @param[in] aMat     Vector.
 * @param[in] WhichValue    Value, which is in the vector.
 *
 * @return  Vector of items found
 */

template<  typename T1, typename T2, typename T3 >
auto
reshape(
        moris::Matrix< T1 , T2 > & aA,
        T3                aB,
        T3                aC)
-> decltype( reshape( aA.matrix_data(),aB,aC) )
{
    return reshape( aA.matrix_data(),aB,aC );
}

}



#endif /* PROJECTS_LINALG_SRC_FN_RESHAPE_HPP_ */
