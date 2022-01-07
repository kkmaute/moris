/*
 * fn_find_unique_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    //Find all indices
    template< typename ET >
    auto
    find_unique( ET & aA )
    ->decltype( arma::find_unique( aA ) )
    {
        return arma::find_unique( aA );
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_FIND_UNIQUE_ARMA_HPP_ */
