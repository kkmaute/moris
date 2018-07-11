#ifndef MORIS_LINALG_FN_DIAG_HPP_
#define MORIS_LINALG_FN_DIAG_HPP_

// MORIS library header files.
#include "fn_isvector.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T>
    moris::Mat< T >
    diag(
        moris::Mat< T > const & aA,
        moris::size_t     const & ak = 0 )
    {
        if (moris::isvector( aA ))
        {
            return arma::diagmat( aA.data() );
        }
        else
        {
            return arma::diagvec( aA.data(), ak );
        }
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T>
    moris::Mat< T >
    diag(
        moris::Mat< T > const & aA,
        moris::size_t     const & ak = 0 )
    {
        if (moris::isvector( aA ))
        {
            return aA.data().asDiagonal();
        }
        else
        {
            return aA.data().diagonal();
        }
    }
}

#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Extract the k-th diagonal from matrix A.
     *
     * @param[in] aA Matrix.
     * @param[in] ak Diagonal index.
     * The argument k is optional; by default the main diagonal
     * is extracted (k=0).\n
     * For k > 0, the k-th super-diagonal is extracted
     * (top-right corner).\n
     * For k < 0, the k-th sub-diagonal is extracted
     * (bottom-left corner).\n
     *
     * @return Creates a vector from the diagonal of a matrix such that
     * @f$ \mathbf{v}_{i}=\mathbf{A}_{ii}@f$ \n
     * The extracted diagonal is interpreted as a column vector.
     */
    template< typename T>
    auto
    diag(
            moris::Mat< T > const & aA,
            moris::size_t     const & ak = 0 )
    ->decltype( moris::Math::diag( aA, ak ) )
    {
        return moris::Math::diag( aA, ak );
    }
}

#endif /* MORIS_FN_DIAG_HPP_ */
