/*
 * fn_sum_Arma.hpp
 *
 *  Created on: Aug 30, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{

    /*template< typename ET, typename Data >
    auto
    sum( const arma::eOp<Data,ET> & aA )
        -> decltype( arma::as_scalar(arma::accu( aA )) )
    {
        return arma::as_scalar(arma::accu( aA ));
    }*/

    template< typename Type>
    Type
    sum( const arma::Mat<Type> & aA )
    {
        return arma::as_scalar(arma::accu( aA ));;
    }

//    template< typename ET >
//    auto
//    sum( ET & aA )
//        -> decltype( arma::as_scalar(arma::accu( aA )) )
//    {
//        return arma::as_scalar(arma::accu( aA ));
//    }
}


#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_SUM_ARMA_HPP_ */
