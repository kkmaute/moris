#ifndef MORIS_LINALG_CL_ARMA_SP_MAT_HPP_
#define MORIS_LINALG_CL_ARMA_SP_MAT_HPP_

// Third-party header files.
#include <armadillo>

// MORIS library header files.
#include "core.hpp"
#include "cl_Base_Arma_Mat.hpp" // LNA/src
#include "cl_Mat.hpp" // LNA/src
#include "fn_isrow.hpp" // LNA/src
#include "fn_trans.hpp" // LNA/src

// Class forward declarations.
namespace moris {
    /**
     * Armadillo sparse matrix wrapper
     */
    template< typename T >
    class Arma_Sp_Mat;
}

template< typename T >
class moris::Arma_Sp_Mat : public moris::Base_Arma_Mat< arma::SpMat< T > >
{

public:

    typedef typename arma::SpMat< T >::elem_type type_t;

    // default constructor
    inline
    Arma_Sp_Mat();

    // -------------------------------------------------------------------------

    // copy constructor
    inline
    Arma_Sp_Mat(
            moris::Arma_Sp_Mat< T > const & mat );

    // -------------------------------------------------------------------------

    // templated constructor
    template< typename A >
    Arma_Sp_Mat(
            A const & X );

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Arma_Sp_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems );

    // -------------------------------------------------------------------------

    // sparse matrix constructor
    inline
    Arma_Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues );

    // -------------------------------------------------------------------------

    // sparse matrix constructor with number of elements
    inline
    Arma_Sp_Mat(
            moris::Mat< moris::uint > const & aRowInd,
            moris::Mat< moris::uint > const & aColInd,
            moris::Mat< T >           const & aValues,
            moris::size_t             const & aIElems,
            moris::size_t             const & aJElems );

    // -------------------------------------------------------------------------

    // default destructor
    inline
    ~Arma_Sp_Mat();

    // -------------------------------------------------------------------------

    // inheriting operator= from parent class
    using moris::Base_Arma_Mat< arma::SpMat< T > >::operator=;

    // -------------------------------------------------------------------------

    // FIXME: This should go in cl_Base_Arma_Mat based on the structure of the
    //        Mat and Sp_Mat classes.

    // two-index access
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index )
    -> decltype( this->mMat( i_index, j_index ) );

    // -------------------------------------------------------------------------

    // const two-index access
    auto
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index ) const
    -> decltype( this->mMat( i_index, j_index ) );

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
     * @brief Dummy call for squeezing unused memory from a sparse matrix
     */
    void
    compress();
};

// -----------------------------------------------------------------------------

// Template implementation file.
#include "cl_Arma_Sp_Mat.tpp"

#endif /* MORIS_LINALG_CL_ARMA_SP_MAT_HPP_ */
