/*
 * fn_join_horiz_Arma.hpp
 *
 *  Created on: Jun 2, 2021
 *      Author: momo
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>


namespace moris
{
    template< typename ET >
    auto
    join_horiz( const ET          & aA,
                const ET          & aB)
    -> decltype( arma::join_horiz( aA, aB ) )
    {

        return arma::join_horiz( aA, aB );

    }


}




#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_FN_JOIN_HORIZ_ARMA_HPP_ */
