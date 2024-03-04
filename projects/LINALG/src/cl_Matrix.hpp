/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Matrix.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_CL_MATRIX_HPP_
#define PROJECTS_LINALG_SRC_CL_MATRIX_HPP_

#include "assert.hpp"
#include "moris_typedefs.hpp"
#include <ostream>

namespace moris
{
    /*
     * This is the Matrix class which all implementations specialize. All
     * functions in this headers throw errors because this is never used.
     * To create a new implementation, look into the implementation directories.
     * In short, these implementations have a specific value of Matrix_Type template
     * which allows the compiler to deduce the correct class.
     *
     * When using moris linear algebra, this header should be included.
     * When adding another implementation, the header should be included at the bottom
     * of this header encapsulated by the correct compiler flags.
     */
    template< typename Matrix_Type >
    class Matrix
    {
        // These member variables are here for error throwing
        Matrix_Type mMatrix;

      public:
        // note for this example this typedef is defined but does nothing
        typedef real Data_Type;

        Matrix()
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // Constructor with no fill value
        Matrix(
                size_t const & aNumRows,
                size_t const & aNumCols )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // Constructor with fill value

        /* Constructor for std::initializer_list<std::initializer_list>,
         * allows for {{1,2,3},{4,5,6},{7,8,9}} to be passed as input
         * and receive,
         * [ 1 2 3 ]
         * [ 4 5 6 ]
         * [ 7 8 9 ]
         *
         * as output
         */
        Matrix( std::initializer_list< std::initializer_list< Data_Type > > const & aInitList )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // template constructor
        template< typename A >
        Matrix( A const & X )
                : mMatrix( X )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // template constructor
        Matrix( Matrix< Matrix_Type > const & X )
                : mMatrix( X )
        {
        }

        // Copy operations
        Matrix< Matrix_Type >
        copy() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix;
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
        void
        resize(
                const size_t& aNumRows,
                const size_t& aNumCols )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
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
        void
        set_size(
                const size_t& aNumRows,
                const size_t& aNumCols )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        void
        set_size(
                const size_t&    aNumRows,
                const size_t&    aNumCols,
                const Data_Type& aFillValue )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Reshape object while preserving elements in the original layout
         *
         * @param[in] aNumRows Number of Rows.
         * @param[in] aNumCols Number of Columns.
         */
        void
        reshape(
                const size_t& aNumRows,
                const size_t& aNumCols )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief In-place / in-situ transpose of matrix X
         */
        void
        inplace_trans()
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
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
        //    void
        //    fill(const Type & aFillValue)
        //    {
        //        MORIS_ERROR(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        //    }

        // -------------------------------------------------------------------------

        void
        set_row(
                size_t                       aRowIndex,
                const Matrix< Matrix_Type >& aRow )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        void
        set_column(
                size_t                       aColumnIndex,
                const Matrix< Matrix_Type >& aColumn )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        void
        get_column(
                size_t                 aColumnIndex,
                Matrix< Matrix_Type >& aColumn ) const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------
        /**
         * Get the number of columns in a data set, similar to Matlab cols().
         *
         * @return Number of columns.
         */
        size_t
        n_cols() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
        }

        /**
         * Get the number of rows in a data set, similar to Matlab rows().
         *
         * @return Number of rows.
         */
        size_t
        n_rows() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
        }

        // -------------------------------------------------------------------------

        /**
         * Returns the length of a vector. Throws error if neither rows nor cols are equal 1.
         */
        size_t
        length() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
        }

        // -------------------------------------------------------------------------

        /**
         * Returns the number of elements in the %matrix.
         *
         * @return Number of elements in the %matrix.
         *
         */
        size_t
        numel() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
        }

        // -------------------------------------------------------------------------

        //    Type*
        //    data() const
        //    {
        //        MORIS_ERROR(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        //        return mDummy;
        //    }

        // -------------------------------------------------------------------------

        Matrix_Type&
        matrix_data()
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix;
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

        Data_Type
        max() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
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

        Data_Type
        max(
                moris::uint& aRowIndex,
                moris::uint& aColIndex )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0.0;
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
        Data_Type
        min() const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0;
        }

        /**
         * Member function.
         * @param[in] aRowIndex Row index of min value of an object.
         * @param[in] aColIndex Column index of min value of an object.
         *
         * @return The extremum value of an object and store the location
         * of the extremum value in the provided variable(s)
         *
         */
        Data_Type
        min(
                moris::uint& aRowIndex,
                moris::uint& aColIndex ) const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return 0.0;
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator()
         *
         * @param[in] aRowIndex Row index for which data should be accessed.
         * @param[in] aColIndex Column index for which data should be accessed.
         */
        Data_Type&
        operator()(
                const size_t& aRowIndex,
                const size_t& aColIndex )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix( aRowIndex, aColIndex );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator()
         *
         * @param[in] aRowIndex Row index for which data should be accessed.
         * @param[in] aColIndex Column index for which data should be accessed.
         */

        const Data_Type&
        operator()(
                const size_t& aRowIndex,
                const size_t& aColIndex ) const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix( aRowIndex, aColIndex );
        }

        // -------------------------------------------------------------------------

        Data_Type&
        operator()( const size_t& aIndex )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix( aIndex );
        }

        // -------------------------------------------------------------------------

        const Data_Type&
        operator()( const size_t& aIndex ) const
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return mMatrix( aIndex );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator =
         *
         * @param[in] moris matrix.
         */

        Matrix< Matrix_Type >&
        operator=( const Matrix< Matrix_Type >& aData )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return *this;
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator =
         *
         * @param[in] Data_Type.
         */

        Matrix< Matrix_Type >&
        operator=( const Data_Type& aData )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return *this;
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator =
         *
         * @param[in] Expression.
         */

        template< typename E >
        Matrix< Matrix_Type >&
        operator=( const E& aExpression )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
            return *this;
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator +=
         *
         * @param[in] moris matrix.
         */

        void
        operator+=( const Matrix< Matrix_Type >& aMatrix )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator +=
         *
         * @param[in] Data_Type.
         */

        void
        operator+=( const Data_Type& aData )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator +=
         *
         * @param[in] Expression.
         */

        template< typename E >
        void
        operator+=( const E& aExpression )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator +=
         *
         * @param[in] moris matrix.
         */

        void
        operator-=( const Matrix< Matrix_Type >& aMatrix )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator -=
         *
         * @param[in] Data_Type.
         */

        void
        operator-=( const Data_Type& aData )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------

        /**
         * @brief Overloaded moris::Matrix_Base::operator -=
         *
         * @param[in] Expression.
         */

        template< typename E >
        void
        operator-=( const E& aExpression )
        {
            MORIS_ERROR( false, "Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?" );
        }

        // -------------------------------------------------------------------------
    };

    template< typename Matrix_Type >
    std::ostream& operator<<( std::ostream& os, const Matrix< Matrix_Type >* aMatrix )
    {
        if ( aMatrix )
        {
            os << *aMatrix;
        }
        else
        {
            os << "NULL";
        }
        return os;
    }

    template< typename Matrix_Type >
    std::ostream& operator<<( std::ostream& os, const Matrix< Matrix_Type >& aMatrix )
    {
        for ( uint i = 0; i < aMatrix.numel(); i++ )
        {
            os << aMatrix( i );
            if ( i < aMatrix.numel() - 1 )
            {
                os << ",";
            }
        }
        return os;
    }
}    // namespace moris

#ifdef MORIS_USE_ARMA
#include "cl_Matrix_Arma_Dynamic.hpp"
#endif

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/cl_Matrix_Eigen_3x3.hpp"
#include "Eigen_Impl/cl_Matrix_Eigen_3x1.hpp"
#include "Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp"
#endif

#include "fn_print.hpp"
#include "op_minus.hpp"
#include "op_plus.hpp"
#include "op_div.hpp"
#include "op_times.hpp"

#endif /* PROJECTS_LINALG_SRC_CL_MATRIX_HPP_ */
