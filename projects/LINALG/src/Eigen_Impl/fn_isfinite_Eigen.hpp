/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_isfinite_Eigen.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_
#include <Eigen/Dense>
#include "fn_Eigen_Aliasing.hpp"

namespace moris
{
    template< typename ET >
    bool isfinite( const Eigen::MatrixBase< ET >& aA )
    {
        if constexpr ( ( Eigen::internal::traits< ET >::Flags & Eigen::DirectAccessBit )
                       && linalg_internal::has_eigen_data< ET >::value )
        {
            const auto& tValues = aA.derived();
            const auto  tCount  = tValues.size();
            if ( tCount == 0 )
            {
                return true;
            }
            if ( tValues.innerStride() == 1
                    && ( tValues.outerSize() == 1 || tValues.outerStride() == tValues.innerSize() ) )
            {
                // allFinite() in Eigen 3.4 tests (A-A) for NaNs, reading each
                // entry repeatedly through expression evaluators at -O0.
                const auto* tData = tValues.data();
                for ( Eigen::Index i = 0; i < tCount; ++i )
                {
                    if ( !(Eigen::numext::isfinite)( tData[ i ] ) )
                    {
                        return false;
                    }
                }
                return true;
            }
        }

        // Preserve native expression nesting and arbitrary strided views.
        return aA.array().isFinite().all();
    }

}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_FN_ISFINITE_EIGEN_HPP_ */
