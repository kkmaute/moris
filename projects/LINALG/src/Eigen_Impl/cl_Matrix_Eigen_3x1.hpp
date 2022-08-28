/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Matrix_Eigen_3x1.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X1_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X1_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template< typename Type >
class Matrix<Eigen::Matrix<Type, 3, 1>>
{
private:
    Eigen::Matrix< Type, 3, 1 > mMatrix;

    void fill_with_NANs()
    {
#ifdef MATRIX_FILL
        if (std::numeric_limits<Data_Type>::has_quiet_NaN)
        {
            mMatrix.fill( std::numeric_limits<Data_Type>::quiet_NaN() );
        }
        else
        {
            mMatrix.fill( std::numeric_limits<Data_Type>::max() );
        }
#endif
    }

public:
    typedef Type Data_Type;

    Matrix()
    {
    };

//    // template constructor
//    template< typename A >
//    Matrix(A const & X ):
//                         mMatrix(X)
//    {
//
//    }

    Matrix( size_t const & aNumRows,
            size_t const & aNumCols)
    {
        MORIS_ASSERT( aNumRows == 3, "Number of rows has to be 3 for 3x1 vector");
        MORIS_ASSERT( aNumCols == 1, "Number of cols has to be 1 for 3x1 vector");

        this->fill_with_NANs();
    }

    Matrix( size_t const & aNumRows,
            size_t const & aNumCols,
            Type const & aFillVal )
    {
        MORIS_ASSERT( aNumRows == 3, "Number of rows has to be 3 for 3x1 vector");
        MORIS_ASSERT( aNumCols == 1, "Number of cols has to be 1 for 3x1 vector");
        mMatrix.fill( aFillVal );
    }

    Matrix( std::initializer_list<std::initializer_list<Type> > const & aInitList )
    {
        size_t i = 0;
        size_t j = 0;

        MORIS_ASSERT( aInitList.size() == 3, "The number of rows have to be 3 for 3x1 vector" );
        MORIS_ASSERT( aInitList.begin()->size() == 1, "The number of columns have to be 1 for 3x1 vector" );

        for(const auto tRow : aInitList) // loop over number of rows
        {
            for(const auto tCol : tRow) // loop over every value in the row
            {
                mMatrix(i, j) = tCol;
                ++j;
            }
            j = 0;
            ++i;
        }
    }

    // template constructor
    Matrix(Eigen::Matrix<Type, 3,1 > const & X ):
        mMatrix(X)
    {

    }

    // template constructor
    template< typename A >
    Matrix(A const & X ):
    mMatrix(X)
    {

    }

    /**
     * Returns a copy of the vector
     *
     * @return Copy of the vector
     */
    Matrix< Eigen::Matrix< Type, 3, 1 > >
    copy() const
    {
        Matrix< Eigen::Matrix< Type, 3, 1 > >tMatCopy;
        tMatCopy.matrix_data() = mMatrix;
        return tMatCopy;
    }

    /**
     * Inserts a given value in all entries of the vector
     *
     * @param[in] aFillValue Given value
     */
    void
    fill( const Type & aFillValue )
    {
        mMatrix.fill( aFillValue );
    }

    /**
     * Get the number of columns = 1 in a data set, similar to Matlab cols().
     *
     * @return Number of columns = 1.
     */
    size_t
    n_cols() const
    {
        return 1;
    }

    /**
     * Get the number of rows = 3 in a data set, similar to Matlab rows().
     *
     * @return Number of rows = 3.
     */
    size_t
    n_rows() const
    {
        return 3;
    }

    /**
     * Returns the number of elements in the %matrix.
     *
     * @return Number of elements in the %matrix.
     *
     */
    size_t
    numel() const
    {
        return mMatrix.size();
    }

//    void set_row( size_t aRowIndex, Matrix< Type, Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic > > const & aRow )
//    {
//        MORIS_ASSERT( aRow.n_rows() == 1, "aRow needs to be a row matrix" );
//        MORIS_ASSERT( aRowIndex < 3, "Specified row index out of bounds" );
//        MORIS_ASSERT( aRow.n_cols() == 1, "Dimension mismatch ( argument matrix has to have have only one column )");
//
//        size_t tROW_INDEX = 0;
//        mMatrix.row(aRowIndex) = aRow.matrix_data().row(tROW_INDEX);
//    }
//
//    void set_column( size_t aColumnIndex,
//                     Matrix<Type, Eigen::Matrix< Type, Eigen::Dynamic, Eigen::Dynamic > > & aColumn )
//    {
//        MORIS_ASSERT( aColumn.n_cols() == 1, "aColumn needs to be a column matrix" );
//        MORIS_ASSERT( aColumnIndex == 0, "Specified column index out of bounds. Has to be 0" );
//        MORIS_ASSERT( aColumn.n_rows() == 3, "Dimension mismatch ( argument matrix has to have 3 rows )");
//
//        size_t tCOLUMN_INDEX = 0;
//        mMatrix.col( aColumnIndex ) = aColumn.matrix_data().col( tCOLUMN_INDEX );
//    }
//
//
//    void
//    get_column( size_t aColumnIndex,
//                Matrix<Type, Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> & aColumn) const
//    {
//        MORIS_ASSERT( aColumn.n_cols() == 1,"aColumn needs to be a column matrix");
//        MORIS_ASSERT( aColumnIndex == 0,"Specified column index out of bounds. Has to be 0");
//        MORIS_ASSERT( aColumn.n_rows() == 3,"Dimension mismatch ( argument matrix has to have 3 rows");
//        const size_t tCOLUMN_INDEX = 0;
//        aColumn.matrix_data().col( tCOLUMN_INDEX ) = mMatrix.col( aColumnIndex );
//    }
//
//    Matrix<Type, Eigen::Matrix<Type, 3, 1 > >
//    get_column( size_t aColumnIndex ) const
//    {
//        MORIS_ASSERT(aColumnIndex == 0,"Specified column index out of bounds. Has to be 0");
//        return mMatrix.col(aColumnIndex);
//    }
//
//    void get_row( size_t aRowIndex, Matrix<Type, Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>> & aRow) const
//    {
//        MORIS_ASSERT( aRow.n_rows() == 1,"aRow needs to be a row matrix" );
//        MORIS_ASSERT( aRowIndex < 3,"Specified row index out of bounds" );
//        MORIS_ASSERT( aRow.n_cols() == 1,"Dimension mismatch ( argument matrix has to have have only one column )" );
//
//        const size_t tROW_INDEX = 0;
//        aRow.mMatrix.row( tROW_INDEX ) = mMatrix.row( aRowIndex );
//    }
//
//    Matrix< Type, Eigen::Matrix< Type, 1, 1 > >
//    get_row(size_t aRowIndex) const
//    {
//        MORIS_ASSERT( aRowIndex < 3,"Specified row index out of bounds" );
//        return mMatrix.row(aRowIndex);
//    }

    const Type*
    data() const
    {
        return mMatrix.data();
    }

    Type*
    data()
    {
        return mMatrix.data();
    }

    inline
    Eigen::Matrix<Type, 3, 1> &
    matrix_data()
    {
        return mMatrix;
    }

    inline
    Eigen::Matrix<Type, 3, 1> const &
    matrix_data() const
    {
        return mMatrix;
    }

    Type
    max() const
    {
        return mMatrix.maxCoeff( );
    }

    Type
    min() const
    {
         return mMatrix.minCoeff( ) ;
    }

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aRowIndex Row index for which data should be accessed.
     * @param[in] aColIndex Column index for which data should be accessed.
     */
    inline
    Type &
    operator()( size_t const & aRowIndex,
                size_t const & aColIndex )
    {
        return mMatrix( aRowIndex, aColIndex );
    }

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aRowIndex Row index for which data should be accessed.
     * @param[in] aColIndex Column index for which data should be accessed.
     */
    const Type &
    operator()(const size_t & aRowIndex,
               const size_t & aColIndex) const
    {
        return mMatrix( aRowIndex, aColIndex );
    }

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aIndex Index for which data should be accessed.
     */
    inline
    Type &
    operator()( const size_t & aIndex )
    {
        return mMatrix( aIndex );
    }

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aIndex Index for which data should be accessed.
     */
    const Type &
    operator()( const size_t & aIndex ) const
    {
        return mMatrix( aIndex );
    }

    /**
     * Returns the length of a vector. Thows error neither rows nor cols are equal 1.
     */
    size_t
    length() const
    {
            return 3;
    }
};
}

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X1_HPP_ */

