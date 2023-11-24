/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INTERSECT_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INTERSECT_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    include "cl_Matrix.hpp"

            namespace moris
    {
        template< typename ET >
        auto
        intersect(
                const ET& aAmat,
                const ET& aBmat )
                ->decltype( arma::intersect( aAmat, aBmat ) )

        {
            return arma::intersect( aAmat, aBmat );
        }

        //-------------------------------------------------------------------------

        template< typename ET, typename ETT >
        void
        intersect(
                const ET& aAmat,
                const ET& aBmat,
                ET&       aCmat,
                ETT&      aAind,
                ETT&      aBind )
        {
            // create temporary vector of uword
            arma::uvec tAind;
            arma::uvec tBind;

            // call arma function
            arma::intersect( aCmat, tAind, tBind, aAmat, aBmat );

            // store temporary vectors in output vectors
            aAind = arma::conv_to< ETT >::from( tAind );
            aBind = arma::conv_to< ETT >::from( tBind );
        }
    }    // namespace moris

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_INTERSECT_EIGEN_HPP_ */
