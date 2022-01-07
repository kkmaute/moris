/*
 * fn_vectorize_Arma.hpp
 *
 *  Created on: May 17, 2021
 *      Author: maute
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_VECTORIZE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_VECTORIZE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename E >
    auto
    vectorize(
            const E & aA,
            bool      aToRowVector=false)
    -> decltype ( arma::vectorise(aA, aToRowVector) )
    {
        return arma::vectorise(aA, aToRowVector);
    }
} /* moris namespace */

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_VECTORIZE_ARMA_HPP_ */
