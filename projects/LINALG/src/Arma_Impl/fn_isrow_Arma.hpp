/*
 * fn_isrow_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    isrow( ET const & aA )
    ->decltype( aA.is_rowvec() )
    {
        return aA.is_rowvec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISROW_ARMA_HPP_ */
