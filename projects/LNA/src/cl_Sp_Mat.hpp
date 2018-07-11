#ifndef MORIS_LINALG_CL_SP_MAT_HPP_
#define MORIS_LINALG_CL_SP_MAT_HPP_

// Third-party header files.
#include <utility>

// MORIS library header files.
#include "typedefs.hpp" // COR/src
#include "cl_Base_Mat.hpp" // LNA/src
#ifdef MORIS_USE_ARMA
#include "cl_Base_Arma_Mat.hpp" // LNA/src
#include "cl_Arma_Sp_Mat.hpp" // LNA/src
#endif

#ifdef MORIS_USE_EIGEN
#include "cl_Base_Eigen_Mat.hpp" // LNA/src
#include "cl_Eigen_Sp_Mat.hpp" // LNA/src
#endif
#include "cl_Mat.hpp" // LNA/src

// Class forward declarations.
namespace moris {
    template< typename T >
    class Sp_Mat;

template< typename T >
#ifdef MORIS_USE_ARMA
    using Sp_Mat_Impl  = Arma_Sp_Mat< T >;
#elif  MORIS_USE_EIGEN
    using Sp_Mat_Impl  = Eigen_Sp_Mat< T >;
#endif
}

/**
 * @brief sparse matrix class
 *
 * moris::Sp_Mat is a template wrapper around arma::SpMat and
 * Eigen::SparseMatrix classes
 *
 * This class converts a sparse or full matrix to sparse form by squeezing out
 * any zero elements
 *
 * Example:
 * @include LNA/src/cl_Sp_Mat/cl_Sp_Mat_form1.inc
 *
 *
 */
template< typename T >
    class moris::Sp_Mat : public moris::Base_Mat< moris::Sp_Mat_Impl< T > >
{

public:

    /**
     * @brief default constructor
     *
     * sparse matrix default constructor
     */
    inline
    Sp_Mat() = default;

    // -------------------------------------------------------------------------

    /**
     * @brief copy constructor
     *
     * sparse matrix copy constructor
     *
     * @param[in] aMat Given constant sparse matrix to be copied
     */
    inline
    Sp_Mat(
            moris::Sp_Mat< T > const & aMat )
        : moris::Base_Mat< moris::Sp_Mat_Impl< T > >( aMat )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * @brief copy constructor
     *
     * sparse matrix copy constructor
     *
     * @param[in] X Given constant sparse matrix to be copied
     */
    template< typename A >
    Sp_Mat(
            A const & X )
        : moris::Base_Mat< moris::Sp_Mat_Impl< T > >( X )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * @brief initialization constructor
     *
     * sparse matrix initialization constructor
     *
     * @param[in] i_elems Number of rows in matrix
     * @param[in] j_elems Number of cols in matrix
     */
    inline
    Sp_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems )
        : moris::Base_Mat< moris::Sp_Mat_Impl< T > >( i_elems, j_elems )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * @brief batch insertion constructor
     *
     * sparse matrix batch insertion constructor
     *
     * @param[in] aRowInd Row indices of matrix
     * @param[in] aColInd Column indices of matrix
     * @param[in] aValues Values of sparse matrix
     *
     * The row and column matrices are dense matrices with size 1xN or Nx1,
     * where N is the number of values to be inserted. The location matrix
     * which is constructed based on row and column matrices is a dense matrix
     * with a size of 2xN or Nx2. The location of the i-th element is specified
     * by the contents of the i-th column or i-th row of the location matrix,
     * where the row is in location (0,i), and the column is in location(1,i) or
     * vise-versa. The values is a dense column vector containing the values to
     * be inserted. It must have the same element type as the sparse matrix. The
     * value in values[i] will be inserted at the location specified by the i-th
     * column of the locations matrix.
     *
     * Example:
     * @include LNA/src/cl_Sp_Mat/cl_Sp_Mat_form1.inc
     */
    inline
    Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues )
        : moris::Base_Mat< Sp_Mat_Impl< T > >(aRowInd, aColInd, aValues)
    {
    }

    // -------------------------------------------------------------------------

    /**
     * @brief batch insertion constructor
     *
     * sparse matrix batch insertion constructor
     *
     * @param[in] aRowInd Row indices of matrix
     * @param[in] aColInd Column indices of matrix
     * @param[in] aValues Values of sparse matrix
     * @param[in] aIElems Number of rows in matrix
     * @param[in] aJElems Number of cols in matrix
     *
     * The row and column matrices are dense matrices with size 1xN or Nx1,
     * where N is the number of values to be inserted. The location matrix
     * which is constructed based on row and column matrices is a dense matrix
     * with a size of 2xN or Nx2. The location of the i-th element is specified
     * by the contents of the i-th column or i-th row of the location matrix,
     * where the row is in location (0,i), and the column is in location(1,i) or
     * vise-versa. The values is a dense column vector containing the values to
     * be inserted. It must have the same element type as the sparse matrix. The
     * value in values[i] will be inserted at the location specified by the i-th
     * column of the locations matrix.
     *
     * Example:
     * @include LNA/src/cl_Sp_Mat/cl_Sp_Mat_form2.inc
     */
    inline
    Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues,
            moris::size_t             const & aIElems,
            moris::size_t             const & aJElems )
        : moris::Base_Mat< moris::Sp_Mat_Impl< T > >( aRowInd, aColInd, aValues, aIElems, aJElems )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * @brief default destructor
     *
     * sparse matrix default destructor.
     */
    inline
    ~Sp_Mat() = default;

    // -------------------------------------------------------------------------

    using moris::Base_Mat< moris::Sp_Mat_Impl< T > >::operator=;

    // -------------------------------------------------------------------------

    /**
     * @brief overload moris sparse matrix operator
     *
     * Overloaded moris::Sp_Mat::operator()
     *
     * @param[in] i_index Row index for which data should be accessed
     * @param[in] j_index Column index for which data should be accessed
     */
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index )
    -> decltype( this->mMat( i_index, j_index ) )
    {
        return this->mMat( i_index, j_index );
    }

    /**
     * @brief overload moris sparse matrix operator for const sparse matrices
     */
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index ) const
    -> decltype( this->mMat( i_index, j_index ) )
    {
        return this->mMat( i_index, j_index );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Changes Matrix Size without preserving data
     *
     * @param[in] aNumRows Number of Rows.
     * @param[in] aNumCols Number of Columns.
     *
     * Changes size of a matrix to aNumRows-by-aNumCols
     * If using Armadillo, the new matrix is always uninitialized.
     * If using Eigen, the new matrix preserves the old values IF the
     * change in size is conservative (e.g. changing a 2-by-3 to a 3-by-2).
     * Otherwise, the new matrix does not preserve the old values.
     */
    auto
    set_size(
        moris::size_t const & aNumRows,
        moris::size_t const & aNumCols )
    -> decltype( this->mMat.set_size( aNumRows, aNumCols ) )
    {
        return this->mMat.set_size( aNumRows, aNumCols );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Returns the number of nonzero elements in a sparse matrix
     *
     * Note : Eigen counts a explicitly assigned 0 entry as a non-zero element.
     *
     * Example:
     * @include LNA/src/cl_Sp_Mat/cl_Sp_Mat_nnz.inc
     */
    moris::uint
    get_nnz()
    {
        return this->mMat.get_nnz();
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Clears the sparsity structure of a sparse matrix
     */
    void
    clear_sparsity()
    {
        this->mMat.clear_sparsity();
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Squeezes unused memory from a sparse matrix.
     *
     * Note : When using armadillo, this is a dummy call because armadillo does
     *        not allow user manipulation of sparse matrix memory.
     */
    void
    compress()
    {
        this->mMat.compress();
    }
};

namespace moris
{
    ///@cond

    // is_Sp_Mat is used template various operators. cond command is used
    // to have doxygen ignore this code, as it should not be part of the API.
    template<class T>
    struct is_Sp_Mat
    {
        static const bool value = false;
    };

    template<class T>
    struct is_Sp_Mat<moris::Sp_Mat<T>>
    {
        static const bool value = true;
    };
    ///@endcond

}

#endif /* MORIS_LINALG_CL_SP_MAT_HPP_ */
