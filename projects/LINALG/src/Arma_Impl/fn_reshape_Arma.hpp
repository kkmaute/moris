/*
 * fn_reshape_Arma.hpp
 *
 *  Created on: Sep 6, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_
#include <armadillo>

namespace moris
{
template<  typename T1, typename T2, typename T3 >
moris::Mat<T1>
reshape(
        moris::Mat< T1 > const & aA,
        T2               const & aB,
        T3               const & aC)
        {
    return arma::conv_to< arma::Mat<T1>>::from(reshape((arma::Mat<T1>)aA.data(),aB,aC));
        }

template<  typename T1, typename T2, typename T3 >
moris::Mat<T1>
reshape(
        moris::Sp_Mat< T1 > const & aA,
        T2               const & aB,
        T3               const & aC)
        {
    return arma::conv_to< arma::Mat<T1>>::from(reshape((arma::Mat<T1>)aA.data(),aB,aC));
        }
}




#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_RESHAPE_ARMA_HPP_ */
