#ifndef MORIS_LINALG_CL_ARMA_MAT_HPP_
#define MORIS_LINALG_CL_ARMA_MAT_HPP_

// Third-party header files.
#include <armadillo>

// MORIS library header files.
#include "core.hpp"
#include "cl_Base_Arma_Mat.hpp" // LNA/src

// Class forward declarations.
namespace moris {

    /**
     * Armadillo matrix wrapper class
     */
    template< typename T >
    class Arma_Mat;
}

template< typename T >
class moris::Arma_Mat : public moris::Base_Arma_Mat< arma::Mat< T > >
{

public:

    // default constructor
    inline
    Arma_Mat();

    // -------------------------------------------------------------------------

    // copy constructor
    inline
    Arma_Mat(
            moris::Arma_Mat< T > const & mat );

    // -------------------------------------------------------------------------

    // template constructor
    template< typename A >
    Arma_Mat(
            A const & X );

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Arma_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems );

    // -------------------------------------------------------------------------

    // array constructor
    inline
    Arma_Mat(
            T*                  & aArray,
            moris::size_t const & aIElems,
            moris::size_t const & aJElems );

    // -------------------------------------------------------------------------

    // initializer constructor
    inline
    Arma_Mat(
            std::initializer_list< std::initializer_list< T > > list);

    // -------------------------------------------------------------------------

    // default destructor
    inline
    ~Arma_Mat();

    // -------------------------------------------------------------------------

    // inhereting operator=
    using moris::Base_Arma_Mat< arma::Mat< T > >::operator=;

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

    // single index access
    auto
    operator()(
            moris::size_t const & i_index )
    -> decltype( this->mMat( i_index ) );

    // -------------------------------------------------------------------------

    // const single index access
    auto
    operator()(
            moris::size_t const & i_index ) const
    -> decltype( this->mMat( i_index ) );

    // -------------------------------------------------------------------------

    // block access
    auto
    operator()(
            std::pair< moris::size_t, moris::size_t > const & aI,
            std::pair< moris::size_t, moris::size_t > const & aJ )
    -> decltype( this->mMat( arma::span( aI.first, aI.second ), arma::span( aJ.first, aJ.second ) ) );

    // -------------------------------------------------------------------------

    // fills all elements
    void
    fill(
            T const & aVal );

    // -------------------------------------------------------------------------

    // fills with an identity matrix
    auto
    eye(
            moris::size_t const & aNumElems )
    -> decltype (this->mMat.eye(aNumElems, aNumElems) );

    // -------------------------------------------------------------------------

    // resizes to given size
    auto
    resize(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( this->mMat.resize( aNumRows, aNumCols ) );

    // -------------------------------------------------------------------------

    // copies size of given object
    auto
    copy_size(
            moris::Arma_Mat< T > const & aMat )
    -> decltype( this->mMat.copy_size( aMat.data() ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    min( ) const
    -> decltype( this->mMat.min( ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    min(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.min( (arma::uword&) aRowIndex, (arma::uword&) aColIndex ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    max( ) const
    -> decltype( this->mMat.max( ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    max(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.max( (arma::uword&) aRowIndex, (arma::uword&) aColIndex ) );

    // -------------------------------------------------------------------------

    // Norm of the matrix.
    moris::real norm();
};

// -----------------------------------------------------------------------------

// Template implementation file.
#include "cl_Arma_Mat.tpp"

#endif /* MORIS_LINALG_CL_ARMA_MAT_HPP_ */
