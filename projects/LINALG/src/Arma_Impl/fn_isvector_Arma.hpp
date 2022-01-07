/*
 * fn_isvector_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isvector( ET const & aA )
    ->decltype( aA.is_vec() )
    {
        return aA.is_vec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISVECTOR_ARMA_HPP_ */
