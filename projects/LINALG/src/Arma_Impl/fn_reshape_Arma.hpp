/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_reshape_Arma.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    // Expose the native lazy operation without adding a competing ADL overload.
    using arma::reshape;
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_ */
