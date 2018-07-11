
#ifndef MORIS_LINALG_OP_GREATER_EQUAL_HPP_
#define MORIS_LINALG_OP_GREATER_EQUAL_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>=(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from(aA.data() >= aB.data());
    }

    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>=(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from( (arma::Mat<T1>) aA.data() >= aB);
    }

    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>=(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from(aA >= (arma::Mat<T2>) aB.data());
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T1, typename T2 >
    auto
    operator>=(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >() )
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        (Eigen::Array< T1   , Eigen::Dynamic, Eigen::Dynamic >) aA.data() >=
        (Eigen::Array< T2   , Eigen::Dynamic, Eigen::Dynamic >) aB.data() );
        return greater_equal.cast< moris::uint >();
    }

    template< typename T1, typename T2 >
    auto
    operator>=(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >())
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        (Eigen::Array< T1   , Eigen::Dynamic, Eigen::Dynamic >) aA.data() >=  aB );
        return greater_equal.cast< moris::uint >();
    }

    template< typename T1, typename T2 >
    auto
    operator>=(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >() )
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        aA >= (Eigen::Array< T2   , Eigen::Dynamic, Eigen::Dynamic >) aB.data() );
        return greater_equal.cast< moris::uint >();
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Determine greater than or equal to
     *
     * @param[in] aA Input matrix
     * @param[in] aB Input matrix
     *
     * This function checks if @f$ \mathbf{A}_{ij} \geq \mathbf{B}_{ij} @f$. It returns a
     * logical matrix with elements set to logical 1 (true) where A is greater than or
     * equal to B; otherwise, it returns logical 0 (false).
     *
     * Example:
     * @include LNA/src/op_greater_equal.inc
     *
     */
    template< typename T1, typename T2 >
    auto
    operator>=(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( moris::Math::operator>=( aA, aB ) )
    {
        return moris::Math::operator>=( aA, aB );
    }

    template< typename T1, typename T2 >
    auto
    operator>=(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    ->decltype( moris::Math::operator>=( aA, aB ) )
    {
        return moris::Math::operator>=( aA, aB );
    }

    template< typename T1, typename T2 >
    auto
    operator>=(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( moris::Math::operator>=( aA, aB ) )
    {
        return moris::Math::operator>=( aA, aB );
    }
}

#endif  /* MORIS_LINALG_OP_GREATER_EQUAL_HPP_ */
