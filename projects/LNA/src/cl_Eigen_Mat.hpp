#ifndef MORIS_LINALG_CL_EIGEN_MAT_HPP_
#define MORIS_LINALG_CL_EIGEN_MAT_HPP_

// Third-party header files.
#include <Eigen>
#define EIGEN_DENSEBASE_PLUGIN "EigenDenseBaseAddons.h"

// MORIS library header files.
#include "core.hpp"
#include "cl_Base_Eigen_Mat.hpp" // LNA/src

// Class forward declarations.
namespace moris {
    /**
     * Eigen Matrix wrapper class
     */
    template< typename T >
    class Eigen_Mat;
}

template< typename T >
class moris::Eigen_Mat : public moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >
{

public:

    // default constructor
    inline
    Eigen_Mat();

    // -------------------------------------------------------------------------

    // copy constructor
    inline
    Eigen_Mat(
            moris::Eigen_Mat< T > const & mat );

    // -------------------------------------------------------------------------

    // template constructor
    template< typename A >
    Eigen_Mat(
            A const & X );

    // -------------------------------------------------------------------------

    // size constructor
    inline
    Eigen_Mat(
            moris::size_t const & i_elems,
            moris::size_t const & j_elems );

    // -------------------------------------------------------------------------

    // array constructor
    inline
    Eigen_Mat(
            T*                  & aArray,
            moris::size_t const & aIElems,
            moris::size_t const & aJElems );

    // -------------------------------------------------------------------------

    // initializer constructor
    inline
    Eigen_Mat(
            std::initializer_list< std::initializer_list< T > > list);

    // -------------------------------------------------------------------------

    // default destructor
    inline
    ~Eigen_Mat();

    // -------------------------------------------------------------------------

    // inhereting operator=
    using moris::Base_Eigen_Mat< Eigen::Matrix< T, Eigen::Dynamic, Eigen::Dynamic > >::operator=;

    // -------------------------------------------------------------------------

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
    -> decltype( this->mMat.block( aI.first, aJ.first, aI.second-aI.first+1,aJ.second-aJ.first+1 ) );

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
    -> decltype (this->mMat.setIdentity(aNumElems, aNumElems) );

    // -------------------------------------------------------------------------

    // resizes to given size
    auto
    resize(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( this->mMat.conservativeResize( aNumRows, aNumCols ) );

    // -------------------------------------------------------------------------

    // copies size of given object
    auto
    copy_size(
            moris::Eigen_Mat< T > const & aMat )
    -> decltype( this->mMat.resizeLike( aMat.data() ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    min( ) const
    -> decltype( this->mMat.minCoeff( ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    min(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.minCoeff( ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    max( ) const
    -> decltype( this->mMat.maxCoeff( ) );

    // -------------------------------------------------------------------------

    // extremum value of an object
    auto
    max(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.maxCoeff( ) );

    // -------------------------------------------------------------------------

    // Norm of the matrix.
    moris::real norm();
};

// -----------------------------------------------------------------------------

// Template implementation file.
#include "cl_Eigen_Mat.tpp"

#endif /* MORIS_LINALG_CL_EIGEN_MAT_HPP_ */
