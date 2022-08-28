/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_comp_abs_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COMP_ABS_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COMP_ABS_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template<typename ET>
    auto comp_abs(ET & aExpTemplate)
    ->decltype(arma::abs(aExpTemplate))
    {
        return arma::abs(aExpTemplate);
    }

    template<typename ET>
    auto comp_abs(ET const & aExpTemplate)
    ->decltype(arma::abs(aExpTemplate))
    {
        return arma::abs(aExpTemplate);
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_COMP_ABS_ARMA_HPP_ */

