#ifndef MORIS_LINALG_FN_GET_SPARSITY_HPP_
#define MORIS_LINALG_FN_GET_SPARSITY_HPP_

// MORIS library header files.
#include "assert.hpp"
#include "cl_Mat.hpp"
#include "cl_Sp_Mat.hpp"

// -----------------------------------------------------------------------------
#ifdef MORIS_USE_ARMA
namespace arma_Math
{
    template< typename T >
    void get_sparsity(
            moris::Mat< moris::uint > & aSparsity,
            const moris::Sp_Mat< T >  & aSpMat )
    {
        auto start = aSpMat.data().begin(); // get iterator to the first element of the sparse matrix
        auto end   = aSpMat.data().end();   // get iterator to the last element of the sparse matrix

        moris::uint ind = 0; // initialize nonzero index

        for( auto it = start; it != end; ++it )
        {
            aSparsity( ind, 0 ) = it.row(); // get the row index for non-zero entry
            aSparsity( ind, 1 ) = it.col(); // get the column index for non-zero entry
            ++ind; // increment nonzero index
        }
    }
}
#endif

// -----------------------------------------------------------------------------

#ifdef MORIS_USE_EIGEN
namespace eigen_Math
{
    template< typename T >
    void get_sparsity(
            moris::Mat< moris::uint > & aSparsity,
            const moris::Sp_Mat< T >  & aSpMat )
    {
        moris::uint ind = 0; // initialize nonzero index

        // loop over non-zero entries column-vise
        for( moris::uint k = 0; k < ( aSpMat.data() ).outerSize(); ++k )
        {
            for( typename Eigen::SparseMatrix< T >::InnerIterator it( aSpMat.data(), k ); it; ++it )
            {
                aSparsity( ind, 0 ) = it.row(); // get the row index for non-zero entry
                aSparsity( ind, 1 ) = it.col(); // get the column index for non-zero entry
                ++ind; // increment nonzero index
            }
        }
    }
}
#endif

// -----------------------------------------------------------------------------

namespace moris
{
    /**
     * @brief Extracts the sparsity structure of a sparse matrix
     *
     * @param[in] aSpMat A sparse matrix.
     *
     * @param[out] aSparsity Sparsity structure of aSpMat. This matrix is of
     *             dimensions n_nonzeros by 2.
     *
     * This function can only be used with sparse matrices. The sparsity is
     * ordered column major.
     *
     * Example:
     * @include "fn_get_sparsity.inc" // snippets LNA
     */
    template< typename T >
    void get_sparsity(
            moris::Mat< moris::uint > & aSparsity,
            const moris::Sp_Mat< T >  & aSpMat )
    {
        // check aSparsity for correct number of columns
        MORIS_ASSERT( aSparsity.n_cols() == 2, "Sparsity structure matrix should be of size n_nnz by 2" );

        moris::Math::get_sparsity( aSparsity, aSpMat );
    }
}

#endif  /* MORIS_LINALG_FN_GET_SPARSITY_HPP_ */
