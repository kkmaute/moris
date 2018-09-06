/*
 * fn_sum_Arma.hpp
 *
 *  Created on: Aug 30, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_

#include "cl_Matrix.hpp"
#include <armadillo>

namespace moris
{

    template< typename ET >
    auto
    sum( const ET & aA )
        -> decltype( arma::accu( aA.matrix_data() ) )
    {
        return arma::accu( aA.matrix_data() );
    }

}


#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_ */
