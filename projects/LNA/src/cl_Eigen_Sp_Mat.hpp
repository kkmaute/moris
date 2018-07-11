#ifndef MORIS_LINALG_CL_EIGEN_SP_MAT_HPP_
#define MORIS_LINALG_CL_EIGEN_SP_MAT_HPP_

// Third-party header files.
#include <Eigen>

// MORIS library header files.
#include "core.hpp"
#include "cl_Base_Eigen_Mat.hpp" // LNA/src
#include "cl_Eigen_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_isrow.hpp" // LNA/src
#include "fn_trans.hpp" // LNA/src

// Class forward declarations.
namespace moris {
    /**
     * Eigen Sparse Matrix wrapper class
     */
    template< typename T >
    class Eigen_Sp_Mat;
}

template< typename T >
class moris::Eigen_Sp_Mat : public moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >
{

public:

    typedef typename Eigen::SparseMatrix< T >::Scalar type_t;

    // default constructor
    inline
    Eigen_Sp_Mat();

    // -------------------------------------------------------------------------

    // copy constructor
    inline
    Eigen_Sp_Mat(
            moris::Eigen_Sp_Mat< T > const & mat );

    // -------------------------------------------------------------------------

    // template constructor
    template< typename A >
    Eigen_Sp_Mat(
            A const & X );

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Eigen_Sp_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems );

    // -------------------------------------------------------------------------

    // sparse matrix constructor
    inline
    Eigen_Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues );

    // -------------------------------------------------------------------------

    // sparse matrix constructor with number of elements
    inline
    Eigen_Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues,
            moris::size_t             const & aIElems,
            moris::size_t             const & aJElems );

    // -------------------------------------------------------------------------

    // default destructor
    inline
    ~Eigen_Sp_Mat();

    // -------------------------------------------------------------------------

    // inheriting operator= from parent class
    using moris::Base_Eigen_Mat< Eigen::SparseMatrix< T > >::operator=;

    // -------------------------------------------------------------------------

    // two-index access
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index )
    -> decltype( this->mMat.coeffRef( i_index, j_index ) );

    // -------------------------------------------------------------------------

    // const two-index access
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index ) const
    -> decltype( this->mMat.coeff( i_index, j_index ) );

    // -------------------------------------------------------------------------

    // get the number of non-zero elements
    moris::uint
    get_nnz();

    // -------------------------------------------------------------------------

    /**
     * @brief Clears the sparsity structure of a sparse matrix
     */
    void
    clear_sparsity();

    // -------------------------------------------------------------------------

    /**
     * @brief Squeezes unused memory from a sparse matrix
     */
    void
    compress();
};

// -----------------------------------------------------------------------------

// Template implementation file.
#include "cl_Eigen_Sp_Mat.tpp"

#endif /* MORIS_LINALG_CL_EIGEN_SP_MAT_HPP_ */
