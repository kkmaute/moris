/*
 * fn_linspace_Eigen.hpp
 *
 *  Created on: Aug 29, 2018
 *      Author: schmidt
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSPACE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSPACE_EIGEN_HPP_
#include <Eigen/Dense>

namespace moris
{
    template< typename T >
    auto
    linspace_base( T            const & aStart,
                   T            const & aEnd,
                   moris::size_t const & aN  )
    -> decltype( Eigen::Matrix< T, Eigen::Dynamic, 1 >::LinSpaced( aN, aStart, aEnd ) )
    {
        typedef Eigen::Matrix< T, Eigen::Dynamic, 1 > tVector;

        return tVector::LinSpaced( aN, aStart, aEnd );
    }
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_LINSPACE_EIGEN_HPP_ */
