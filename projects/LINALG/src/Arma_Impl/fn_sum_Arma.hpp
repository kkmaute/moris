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

    template< typename T >
    auto
    sum( T & aA )
    -> decltype( arma::accu(aA) )
    {
        std::cout<<arma::accu(aA)<<std::endl;

        return arma::accu(aA);
    }

}


#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_ */
