/*
 * fn_dot_Arma.hpp
 *
 *  Created on: Aug 31, 2018
 *      Author: sonne
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DOT_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DOT_ARMA_HPP_

#include "cl_Matrix.hpp"

namespace moris
{
    template< typename ET >
    auto
    dot(
            const ET & aA,
            const ET & aB)
    -> decltype( arma::dot(aA , aB ) )
    {
        return arma::dot( aA , aB );
    }

    template< typename Matrix_Type, typename ET >
    auto
    dot(
            Matrix< Matrix_Type > const & aA,
            ET                    const & aB)
    ->decltype(arma::dot(aA.matrix_data(), aB))
    {
        return arma::dot(aA.matrix_data(), aB);
    }

    template<typename Matrix_Type, typename ET >
    auto
    dot(
            ET                    const & aA,
            Matrix< Matrix_Type > const & aB)
    ->decltype(arma::dot(aA, aB.matrix_data()))
    {
        return arma::dot(aA, aB.matrix_data());
    }
}

#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_DOT_ARMA_HPP_ */
