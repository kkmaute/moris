/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linspace.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_LINSPACE_HPP_
#define PROJECTS_LINALG_SRC_FN_LINSPACE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_linspace_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_linspace_Arma.hpp"
#endif

namespace moris
{
/**
 * @brief Generate a vector with N elements; the values of the elements
 * linearly increase from start to (and including) end.
 *
 *@param[in] aStart Start of vector
 *@param[in] aEnd End pf vector
 *@param[in] aN Number of elements in vector
 *
 */
    template< typename T >
    auto
    linspace( T             const & aStart,
              T             const & aEnd,
              moris::size_t const & aN )
    -> decltype( linspace_base( aStart, aEnd, aN ) )
    {
        return linspace_base( aStart, aEnd, aN );
    }
}

#endif /* PROJECTS_LINALG_SRC_FN_LINSPACE_HPP_ */

