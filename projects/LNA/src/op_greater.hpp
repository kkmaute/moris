#ifndef MORIS_LINALG_OP_GREATER_HPP_
#define MORIS_LINALG_OP_GREATER_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>(
             moris::Mat< T1 > const & aA,
             moris::Mat< T2 > const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from(aA.data() > aB.data());
    }

    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from( (arma::Mat<T1>) aA.data() > aB);
    }

    template< typename T1, typename T2 >
    moris::Mat< moris::uint >
    operator>(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    {
        return arma::conv_to< arma::Mat< moris::uint > >::from(aA > (arma::Mat<T2>) aB.data());
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T1, typename T2 >
    auto
    operator>(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >() )
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        (Eigen::Array< T1   , Eigen::Dynamic, Eigen::Dynamic >) aA.data() >
        (Eigen::Array< T2   , Eigen::Dynamic, Eigen::Dynamic >) aB.data() );
        return greater_equal.cast< moris::uint >();
    }

    template< typename T1, typename T2 >
    auto
    operator>(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >())
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        (Eigen::Array< T1   , Eigen::Dynamic, Eigen::Dynamic >) aA.data() > aB );
        return greater_equal.cast< moris::uint >();
    }

    template< typename T1, typename T2 >
    auto
    operator>(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( Eigen::Matrix< moris::uint, Eigen::Dynamic, Eigen::Dynamic >() )
    {
        Eigen::Matrix< bool, Eigen::Dynamic, Eigen::Dynamic > greater_equal(
        aA > (Eigen::Array< T2   , Eigen::Dynamic, Eigen::Dynamic >) aB.data() );
        return greater_equal.cast< moris::uint >();
    }
}
#endif
// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Element-wise check if one matrix is greater than another.
     *
     * @param[in] aA Matrix.
     * @param[in] aB Matrix.
     *
     * @return Checks if: @f$ \mathbf{A}_{ij} >  \mathbf{B}_{ij} @f$ \n
     *         Element-wise equality evaluation of two objects; generates a
     *         matrix with entries that indicate whether at a given position
     *         the two elements from the two objects are greater (1) or not (0).
    */
    template< typename T1, typename T2 >
    auto
    operator>(
            moris::Mat< T1 > const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( moris::Math::operator>( aA, aB ) )
    {
        return moris::Math::operator>( aA, aB );
    }

    template< typename T1, typename T2 >
    auto
    operator>(
            moris::Mat< T1 > const & aA,
            T2               const & aB )
    ->decltype( moris::Math::operator>( aA, aB ) )
    {
        return moris::Math::operator>( aA, aB );
    }

    template< typename T1, typename T2 >
    auto
    operator>(
            T1               const & aA,
            moris::Mat< T2 > const & aB )
    ->decltype( moris::Math::operator>( aA, aB ) )
    {
        return moris::Math::operator>( aA, aB );
    }
}

#endif  /* MORIS_LINALG_OP_GREATER_HPP_ */
