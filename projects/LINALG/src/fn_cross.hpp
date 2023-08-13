/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_cross.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_CROSS_HPP_
#define PROJECTS_LINALG_SRC_FN_CROSS_HPP_

#include "cl_Matrix.hpp"
#include "fn_isvector.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_cross_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/fn_cross_Arma.hpp"
#endif

namespace moris
{

/**
 * Take the cross product of vectors aA with vector aB
 *  aA x aB. This works for strictly vectors only with Eigen. Arma
 *  seems a little more flexible
 *
 *  @param[in] aA - Vector 1
 *  @param[in] aB - Vector 2
 *  @param[out] aA x aB
 */
template< typename Matrix_Type >
auto
cross(Matrix< Matrix_Type > const & aA,
      Matrix< Matrix_Type > const & aB)
-> decltype( cross( aA.matrix_data(),aB.matrix_data()) )
{
    MORIS_ASSERT(isvector(aA),"Provided matrix aA in cross product is not a vector");
    MORIS_ASSERT(isvector(aB),"Provided matrix aB in cross product is not a vector");
    return cross( aA.matrix_data(), aB.matrix_data()) ;
}
}

#endif /* PROJECTS_LINALG_SRC_FN_CROSS_HPP_ */

