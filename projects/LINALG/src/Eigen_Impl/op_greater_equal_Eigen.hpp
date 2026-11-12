/* Copyright (c) 2026 University of Colorado. Licensed under the MIT license. */
#pragma once

#include "cl_Matrix.hpp"
#include <Eigen/Dense>

namespace Eigen
{
    template< class A, class B >
    auto operator>=( const Eigen::MatrixBase< A >& a, const Eigen::MatrixBase< B >& b ) -> decltype( a.array() >= b.array() )
    {
        return a.array() >= b.array();
    }

    template< class E >
    auto operator>=( const Eigen::MatrixBase< E >& a, const typename E::Scalar& b ) -> decltype( a.array() >= b )
    {
        return a.array() >= b;
    }

    template< class E >
    auto operator>=( const typename E::Scalar& a, const Eigen::MatrixBase< E >& b ) -> decltype( a >= b.array() )
    {
        return a >= b.array();
    }
}    // namespace Eigen
namespace moris
{
    using Eigen::operator>=;
}
