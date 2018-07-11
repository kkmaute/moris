#ifndef MORIS_LINALG_FN_LINSPACE_HPP_
#define MORIS_LINALG_FN_LINSPACE_HPP_

// MORIS library header files.
#include "cl_Mat.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T>
    auto
    linspace(
            T           const & aStart,
            T           const & aEnd,
            moris::size_t const & aN )
    -> decltype( arma::linspace< arma::Mat< T > >( aStart, aEnd, aN ) )
    {
        return arma::linspace< arma::Mat< T > >( aStart, aEnd, aN );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T>
    auto
    linspace(
            T           const & aStart,
            T           const & aEnd,
            moris::size_t const & aN )
    -> decltype( Eigen::Matrix< T, Eigen::Dynamic, 1 >::LinSpaced( aN, aStart, aEnd ) )
    {
        typedef Eigen::Matrix< T, Eigen::Dynamic, 1 > tVector;

        return tVector::LinSpaced( aN, aStart, aEnd ); // eigen doesn't work for unsigned long int.
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Generate a vector with N elements; the values of the elements
     * linearly increase from start to (and including) end.
     *
     *@param[in] aStart Start of vector
     *@param[in] aEnd End pf vector
     *@param[in] aN Number of elements in vector
     *
     * Example:
     * @include LNA/src/fn_linspace.inc
     *
     */
    template< typename T>
//    template< typename T, bool M = ( std::is_same< T, moris::uint >::value ||
//            std::is_same< T, moris::cplx >::value), typename std::enable_if< !M >::type* = nullptr >
    auto
    linspace(
            T           const & aStart,
            T           const & aEnd,
            moris::size_t const & aN )
    -> decltype( moris::Math::linspace( aStart, aEnd, aN ) )
    {
        return moris::Math::linspace( aStart, aEnd, aN );
    }
}

#endif  /* MORIS_LINALG_FN_LINSPACE_HPP_ */
