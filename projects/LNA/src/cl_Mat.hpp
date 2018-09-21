#ifndef MORIS_LINALG_CL_MAT_HPP_
#define MORIS_LINALG_CL_MAT_HPP_

// Third-party header files.
#include <utility>
#include <iostream>

// MORIS library header files.
#include "typedefs.hpp"
#include "cl_Base_Mat.hpp"

#ifdef MORIS_USE_ARMA
#include "cl_Base_Arma_Mat.hpp"
#include "cl_Arma_Mat.hpp"
#endif

#ifdef MORIS_USE_EIGEN
#include "cl_Base_Eigen_Mat.hpp"
#include "cl_Eigen_Mat.hpp"
#endif

#include "op_ostream.hpp"

// Class forward declarations.
namespace moris
{
    template< typename T >
    class Mat;

    template< typename T >
#ifdef MORIS_USE_ARMA
    using Mat_Impl = Arma_Mat< T >;
#elif  MORIS_USE_EIGEN
    using Mat_Impl = Eigen_Mat< T >;
#endif
}

/**
 * The matrix class.
 *
 * A matrix can be declared as follows:
 *
 * @include cl_Mat/Mat.inc
 */
template< typename T >
    class moris::Matrix< DDRMat > : public moris::Base_Mat< moris::Mat_Impl< T > >
{

public:

    /**
     * Mat default constructor.
     */
    inline
    Mat() = default;

    // -------------------------------------------------------------------------

    /**
     * Mat copy constructor.
     */
    inline
    Mat(
            moris::Mat< T > const & mat )
        : moris::Base_Mat< moris::Mat_Impl< T > >( mat )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Mat copy constructor.
     */
    template< typename A >
    Mat(
            A const & X )
        : moris::Base_Mat< moris::Mat_Impl< T > >( X )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Mat initialization constructor.
     *
     * @param[in] aIElems Number of rows in matrix.
     * @param[in] aJElems Number of cols in matrix.
     */
    inline
    Mat(
            moris::size_t const & aIElems,
            moris::size_t const & aJElems )
        : moris::Base_Mat< moris::Mat_Impl< T > >(aIElems, aJElems)
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Mat initialization constructor.
     *
     * @param[in] aArray  Pre-defined array(pointer).
     * @param[in] aIElems Number of rows in matrix.
     * @param[in] aJElems Number of cols in matrix.
     *
     * @note The matrix correspond to a column wise arrangement of the
     * pointer.
     */
    inline
    Mat(
            T*                  & aArray,
            moris::size_t const & aIElems,
            moris::size_t const & aJElems )
        : moris::Base_Mat< moris::Mat_Impl< T > >(aIElems, aJElems)
    {
        this->mMat = moris::Mat_Impl<T>(aArray, aIElems, aJElems);
    }

    // -------------------------------------------------------------------------

    /**
     * Mat initialization constructor.
     *
     * @param[in] aIElems Number of rows in matrix.
     * @param[in] aJElems Number of cols in matrix.
     * @param[in] aVal A given value to set all elements to that.
     *
     */
    inline
    Mat(
            moris::size_t const & aIElems,
            moris::size_t const & aJElems,
            T             const & aVal )
        : moris::Base_Mat< moris::Mat_Impl< T > >(aIElems, aJElems)
    {
        this->mMat.fill( aVal );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Mat initialization constructor
     *
     * Will not work (and doesn't make sense) for empty matrices
     *
     * @note a 2D initialization list is used to construct the matrices.
     * For example - { {1,2,3}, {4,5,6} } will give a 2 by 3 matrix.
     * list.size() corresponds to the number of rows. list.begin()->size()
     * corresponds to the number of columns. We use 'range based for loops'
     * to ensure that the constructor works well for complex numbers too.
     */
    inline
    Mat(
            std::initializer_list< std::initializer_list< T > > list)
        : moris::Base_Mat< moris::Mat_Impl< T > >( list )
    {
    }

    // -------------------------------------------------------------------------

    /**
     * Mat destructor.
     */
    inline
    ~Mat() = default;

    // -------------------------------------------------------------------------

    using moris::Base_Mat< moris::Mat_Impl< T > >::operator=;

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Mat::operator()
     *
     * @param[in] i_index Row index for which data should be accessed.
     * @param[in] j_index Column index for which data should be accessed.
     */
    T &
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index )
//    -> decltype( this->mMat( i_index, j_index ) )
    {
        return this->mMat( i_index, j_index );
    }

    /**
     * @brief Overloaded moris::Mat::operator() for constant moris::Mat
     */
    T const &
    operator()(
            moris::size_t const & i_index,
            moris::size_t const & j_index ) const
//    -> decltype( this->mMat( i_index, j_index ) )
    {
        return this->mMat( i_index, j_index );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded single-index moris::Mat::operator()
     *
     * Grants 0-based column-oriented (i.e. column-by-column) single index data access.
     *
     * @param[in] i_index 0-based index
     *
     */
    auto
    operator()(
            moris::size_t const & i_index )
    -> decltype( this->mMat( i_index ) )
    {
        return this->mMat( i_index );
    }

    // -------------------------------------------------------------------------

    auto
    operator()(
            moris::size_t const & i_index ) const
    -> decltype( this->mMat( i_index ) )
    {
        return this->mMat( i_index );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Mat::operator().
     *
     * @param[in] aI Span of rows for which data should be accessed.
     * @param[in] aJ Span of columns for which data should be accessed.
     *
     * Example:
     * @include cl_Mat/Mat_span.inc
     */
    auto
    operator()(
            std::pair< moris::size_t, moris::size_t > const & aI,
            std::pair< moris::size_t, moris::size_t > const & aJ )
    -> decltype( this->mMat( aI, aJ ) )
    {
        return this->mMat( aI, aJ );
    }

    // -------------------------------------------------------------------------

    /**
     * Initialization operator.
     *
     * @param[in] aVal A given value to set all elements to that.
     *
     * Example:
     * @include cl_Mat/Mat_fill.inc
     */
    void
    fill(
            T const & aVal )
    {
        this->mMat.fill( aVal );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Identity Matrix
     *
     * Converts this matrix into an identity matrix.
     *
     * @param[in] aNumElems Desired size of identity matrix.
     *
     * @return Reference to this matrix.
     */
    auto
    eye(
                moris::size_t const & aNumElems )
    -> decltype( this->mMat.eye(aNumElems) )
    {
        return this->mMat.eye(aNumElems);
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Resizes object while preserving elements in the original layout
     *
     * @param[in] aNumRows Number of Rows.
     * @param[in] aNumCols Number of Columns.
     *
     * Changes size of a matrix to to aNumRows-by-aNumCols
     * while preserving the elements in original layout. If you wish
     * to resize the matrix and do not care about the elements, use moris::set_size,
     * which is faster than resize.
     */
    auto
    resize(
            moris::size_t const & aNumRows,
            moris::size_t const & aNumCols )
    -> decltype( this->mMat.resize( aNumRows, aNumCols ) )
    {
        return this->mMat.resize( aNumRows, aNumCols );
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Resizes object to match that of the input object
     *
     * @param[in] aMat Matrix whose size will be copied
     *
     * Changes size of a matrix to match that the input matrix.
     */
    auto
    copy_size(
            moris::Mat< T > const & aMat )
    -> decltype( this->mMat.copy_size( aMat.data() ) )
    {
        {
            return this->mMat.copy_size( aMat.data() );
        }
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
     * @brief Changes Matrix Size without preserving data
     *
     * @param[in] aNumRows Number of Rows.
     * @param[in] aNumCols Number of Columns.
     * @param[in] aVal A given value to set all elements to that.
     *
     * Changes size of a matrix to aNumRows-by-aNumCols and fills it
     * with aVal.
     */
    auto
    set_size(
        moris::size_t const & aNumRows,
        moris::size_t const & aNumCols,
        T             const & aVal)
    -> decltype( this->mMat )
    {
        this->mMat.set_size( aNumRows, aNumCols );
        this->mMat.fill( aVal );
        return this->mMat;
    }

    // -------------------------------------------------------------------------

    /**
     * Member function.
     *
     * @return The extremum value of an object.
     *
     * Example:
     * @include cl_Mat/Mat_min.inc
     */
    auto
    min() const
    -> decltype( this->mMat.min( ) )
    {
        return this->mMat.min( );
    }

    // -------------------------------------------------------------------------

    /**
     * Member function.
     * @param[in] aRowIndex Row index of min value of an object.
     * @param[in] aColIndex Column index of min value of an object.
     *
     * @return The extremum value of an object and store the location
     * of the extremum value in the provided variable(s)
     *
     */
    auto
    min(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.min( aRowIndex, aColIndex ) )
    {
        return this->mMat.min( aRowIndex, aColIndex );
    }

    // -------------------------------------------------------------------------

    /**
     * Returns minimum value and its location as a tuple.
     *
     * @return The minimum value of this matrix and its location
     * (row index, column index), in that order, stored in a
     * moris::Tuple.
     *
     * @usesTuple
     *
     */
    moris::Tuple< T, moris::uint, moris::uint >
    mint()
    {
        moris::uint rowInd;
        moris::uint colInd;
        T           minVal = this->min(rowInd,colInd);

        return moris::make_tuple(minVal,rowInd,colInd);
    }

    // -------------------------------------------------------------------------

    /**
     * Member function.
     *
     * @return The extremum value of an object.
     *
     * Example:
     * @include cl_Mat/Mat_max.inc
     */
    auto
    max() const
    -> decltype( this->mMat.max( ) ) const
    {
        return this->mMat.max( );
    }

    // -------------------------------------------------------------------------

    /**
     * @param[in] aRowIndex Row index of max value of an object.
     * @param[in] aColIndex Column index of max value of an object.
     *
     * @return The extremum value of an object and store the location
     * of the extremum value in the provided variable(s)
     *
     */
    auto
    max(
            moris::uint & aRowIndex,
            moris::uint & aColIndex ) const
    -> decltype( this->mMat.max( aRowIndex, aColIndex ) )
    {
        return this->mMat.max(  aRowIndex, aColIndex );
    }

    // -------------------------------------------------------------------------

    /**
     * Returns maximum value and its location as a tuple.
     *
     * @return The maximum value of this matrix and its location
     * (row index, column index), in that order, stored in a
     * moris::Tuple.
     *
     * @usesTuple
     *
     */
    moris::Tuple< T, moris::uint, moris::uint >
    maxt()
    {
        moris::uint rowInd;
        moris::uint colInd;
        T           maxVal = this->max(rowInd,colInd);

        return moris::make_tuple(maxVal,rowInd,colInd);
    }

    // -------------------------------------------------------------------------

};

namespace moris
{
    ///@cond

    // is_Mat is used template various operators. cond command is used
    // to have doxygen ignore this code, as it should not be part of the API.
    template<class T>
    struct is_Mat
    {
        static const bool value = false;
    };

    template<class T>
    struct is_Mat<moris::Mat<T>>
    {
        static const bool value = true;
    };
    ///@endcond
}

#endif /* MORIS_LINALG_CL_MAT_HPP_ */
