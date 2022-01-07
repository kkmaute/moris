/*
 * fn_iscol_Arma.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    iscol( ET const & aA )
    ->decltype( aA.is_colvec() )
    {
        return aA.is_colvec();
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_ISCOL_ARMA_HPP_ */
