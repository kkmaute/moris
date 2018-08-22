#ifndef MORIS_LINALG_OP_ELEMWISE_DIV_HPP_
#define MORIS_LINALG_OP_ELEMWISE_DIV_HPP_

// Third-party header files.
// MORIS library header files.
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"

#ifdef MORIS_USE_ARMA
#include <armadillo>
namespace arma_Math
{
    template< typename T1, typename T2 >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            moris::Mat< T2 > const & B )
    -> decltype( A.data() / B.data() )
    {
        return A.data() / B.data();
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            T2               const & B )
    -> decltype( A.data() / B )
    {
        return A.data() / B;
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            T1               const & A,
            moris::Mat< T2 > const & B )
    -> decltype( A / B.data() )
    {
        return A / B.data();
    }
}

#endif

// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
#include <Eigen/Core>

namespace eigen_Math
{
    template< typename T1, typename T2 >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            moris::Mat< T2 > const & B )
    -> decltype( A.data().cwiseQuotient( B.data() ) )
    {
        return A.data().cwiseQuotient( B.data() );
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            T2               const & B )
    -> decltype( A.data() / B )
    {
        return A.data() / B;
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            T1               const & A,
            moris::Mat< T2 > const & B )
    -> decltype( A / B.data() )
    {
        return A / B.data();
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Element wise division operator.
     *
     * @param[in] A Elements of A are the dividend.
     * @param[in] B Elements of B are the divisor.
     *
     * @return Creates a matrix corresponding to element-wise
     * division of the two input matrices.
     *USE_EIGEN
     * Example:
     * @include LNA/src/op_elemwise_div.inc
     *
     */
    template< typename T1, typename T2 >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            moris::Mat< T2 > const & B )
    -> decltype( moris::Math::operator/( A, B ) )
    {
        return moris::Math::operator/( A, B );
    }

    template< typename T1, typename T2 >
    auto
    operator/(
            moris::Sp_Mat< T1 > const & A,
            moris::Sp_Mat< T2 > const & B )
    -> decltype( moris::Math::operator/( A, B ) ) = delete;

    /**
     * @brief Element wise division operator.
     *
     * @param[in] A Elements of A are the dividend.
     * @param[in] B Elements of B are the divisor.
     *
     * @return Creates a matrix corresponding to element-wise
     * division of the matrix A by constant B.
     */
    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            moris::Mat< T1 > const & A,
            T2               const & B )
    -> decltype( moris::Math::operator/( A, B ) )
    {
        return moris::Math::operator/( A, B );
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Sp_Mat< T2 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            moris::Sp_Mat< T1 > const & A,
            T2                  const & B )
    -> decltype( moris::Math::operator/( A, B ) ) = delete;

    /**
     * @brief Element wise division operator.
     *
     * @param[in] A Elements of A are the dividend.
     * @param[in] B Elements of B are the divisor.
     *
     * @return Creates a matrix corresponding to element-wise
     * division of the matrix A by constant B.
     */
    template< typename T1, typename T2,
    bool M = ( moris::is_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            T1               const & A,
            moris::Mat< T2 > const & B )
    -> decltype( moris::Math::operator/( A, B ) )
    {
        return moris::Math::operator/( A, B );
    }

    template< typename T1, typename T2,
    bool M = ( moris::is_Sp_Mat< T1 >::value ),
    typename std::enable_if< ! M >::type* = nullptr >
    auto
    operator/(
            T1                  const & A,
            moris::Sp_Mat< T2 > const & B )
    -> decltype( moris::Math::operator/( A, B ) ) = delete;
}

#endif  /* MORIS_LINALG_OP_ELEMWISE_DIV_HPP_ */
