/*
 * fn_rank_Arma.hpp
 *
 *  Created on: May 30, 2021
 *      Author: momo
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_

#include "cl_Matrix.hpp"
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>


namespace moris
{
    template< typename ET>
    auto
    rank(ET const & aA)
    -> decltype( arma::rank( aA ) )
    {
        return arma::rank( aA );
    }

}



#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RANK_ARMA_HPP_ */
