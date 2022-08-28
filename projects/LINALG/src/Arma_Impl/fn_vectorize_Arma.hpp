/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_vectorize_Arma.hpp
 *
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

