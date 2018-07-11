#ifndef MORIS_LINALG_FN_EYE_HPP_
#define MORIS_LINALG_FN_EYE_HPP_

// MORIS library header files.
#include "cl_Mat.hpp" // LNA/src
#include "core.hpp"

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    inline
    auto
    eye(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( arma::eye( aNumRows, aNumCols ) )
    {
        return arma::eye( aNumRows, aNumCols );
    }
}
#endif

// ----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    inline
    auto
    eye(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( Eigen::MatrixXd::Identity( aNumRows, aNumCols ) )
    {
        return Eigen::MatrixXd::Identity( aNumRows, aNumCols );
    }
}
#endif

// ----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Creates an Identity Matrix
     *
     * @param[in] aNumRows Number of Rows.
     * @param[in] aNumCols Number of Columns.
     * An identity matrix is generated when aNumRows = aNumCols
     *
     * @return Creates and identity matrix @f$ \mathbf{I}_{n}@f$
     * such that @f$ \mathbf{({I}_{n})}_{ij} =  \mathbf{\delta}_{ij} @f$\n
     * Generates a matrix with the elements along the main diagonal
     * set to one and off-diagonal elements set to zero.
     */
    inline
    auto
    eye(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( moris::Math::eye( aNumRows, aNumCols ) )
    {
        return moris::Math::eye( aNumRows, aNumCols );
    }
}

#endif /* MORIS_LINALG_FN_EYE_HPP_ */
