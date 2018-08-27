/*
 * cl_Matrix.hpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_CL_MATRIX_HPP_
#define PROJECTS_LINALG_SRC_CL_MATRIX_HPP_

#include <memory>
#include "assert.hpp"
#include "typedefs.hpp"

namespace moris
{

template<typename Type, typename Matrix_Type>
class Mat_New
{
    Matrix_Type mMatrix;

public:
    Mat_New(){}

    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols)
    {

    }

    // template constructor
    template< typename A >
    Mat_New(A const & X ):
                mMatrix(X)
     {

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
//    void
//    resize(const size_t & aNumRows,
//           const size_t & aNumCols)
//    {
//        mBaseMat->resize(aNumRows,aNumCols);
//    }

//    // -------------------------------------------------------------------------
//
//    /**
//     * @brief Changes Matrix Size without preserving data
//     *
//     * @param[in] aNumRows Number of Rows.
//     * @param[in] aNumCols Number of Columns.
//     *
//     * Changes size of a matrix to aNumRows-by-aNumCols
//     * If using Armadillo, the new matrix is always uninitialized.
//     * If using Eigen, the new matrix preserves the old values IF the
//     * change in size is conservative (e.g. changing a 2-by-3 to a 3-by-2).
//     * Otherwise, the new matrix does not preserve the old values.
//     */
//    void set_size( size_t aNumRows,
//                   size_t aNumColumns )
//    {
//        mBaseMat->set_size(aNumRows,aNumColumns);
//    }
//
//    // -------------------------------------------------------------------------
//
//    /**
//     * Initialization operator.
//     *
//     * @param[in] aVal A given value to set all elements to that.
//     *
//     * Example:
//     * @include cl_Mat/Mat_fill.inc
//     */
//    void
//    fill(const Type & aFillValue)
//    {
//        mBaseMat->fill(aFillValue);
//    }
//
//    // -------------------------------------------------------------------------
//
//    void
//    set_row( size_t aRowIndex,
//             const moris::Matrix_Base<Type, Matrix_Type> & aRow )
//    {
//        mBaseMat->set_row(aRowIndex,aRow);
//    }
//
//    // -------------------------------------------------------------------------
//
//    void
//    set_column( size_t aColumnIndex,
//                const moris::Matrix_Base<Type, Matrix_Type> & aColumn )
//    {
//        mBaseMat->set_column(aColumnIndex,aColumn);
//    }
//
//    // -------------------------------------------------------------------------
//
//    void
//    get_column(size_t aColumnIndex,
//               Matrix_Base<Type, Matrix_Type> & aColumn) const
//    {
//        mBaseMat->get_columns(aColumnIndex,aColumn);
//    }
//
//    // -------------------------------------------------------------------------
//    /**
//     * Get the number of columns in a data set, similar to Matlab cols().
//     *
//     * @return Number of columns.
//     */
//    size_t
//    n_cols()
//    {
//        return mBaseMat->n_cols();
//    }
//
//
//    // -------------------------------------------------------------------------
//
//    /**
//     * Returns the number of elements in the %matrix.
//     *
//     * @return Number of elements in the %matrix.
//     *
//     */
//    size_t
//    numel() const
//    {
//        return mBaseMat->numel();
//    }
//
//    // -------------------------------------------------------------------------
//
//    Type*
//    data() const
//    {
//        return mBaseMat->data();
//    }
//
//    // -------------------------------------------------------------------------
//
    Matrix_Type &
    matrix_data()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented?");
        return mMatrix;
    }
//
//
//    // -------------------------------------------------------------------------
//
//    /**
//     * Member function.
//     *
//     * @return The extremum value of an object.
//     *
//     * Example:
//     * @include cl_Mat/Mat_max.inc
//     */
//
//    Type
//    max() const
//    {
//        return mBaseMat->max();
//    }
//
//    // -------------------------------------------------------------------------
//
//    /**
//     * @param[in] aRowIndex Row index of max value of an object.
//     * @param[in] aColIndex Column index of max value of an object.
//     *
//     * @return The extremum value of an object and store the location
//     * of the extremum value in the provided variable(s)
//     *
//     */
//
//
//    Type
//    max( moris::uint & aRowIndex,
//         moris::uint & aColIndex )
//    {
//        return mBaseMat->max(aRowIndex,aColIndex);
//    }
//
//    // -------------------------------------------------------------------------
//
//    /**
//     * Member function.
//     *
//     * @return The extremum value of an object.
//     *
//     * Example:
//     * @include cl_Mat/Mat_min.inc
//     */
//    Type
//    min() const
//    {
//        return mBaseMat->min();
//    }
//
//    /**
//     * Member function.
//     * @param[in] aRowIndex Row index of min value of an object.
//     * @param[in] aColIndex Column index of min value of an object.
//     *
//     * @return The extremum value of an object and store the location
//     * of the extremum value in the provided variable(s)
//     *
//     */
//    Type
//    min( moris::uint & aRowIndex,
//         moris::uint & aColIndex ) const
//    {
//        return mBaseMat->min(aRowIndex,aColIndex);
//    }
//
    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aRowIndex Row index for which data should be accessed.
     * @param[in] aColIndex Column index for which data should be accessed.
     */
    Type &
    operator()(const size_t & aRowIndex,
               const size_t & aColIndex)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented?");
        return mMatrix(aRowIndex,aColIndex);

    }

    // -------------------------------------------------------------------------

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
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented?");
        return mMatrix(aRowIndex,aColIndex);
    }

};
}

#ifdef MORIS_USE_ARMA
#include "Arma_Impl/cl_Matrix_Arma_Dynamic.hpp"
#endif

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/cl_Matrix_Eigen_3x3.hpp"
#include "Eigen_Impl/cl_Matrix_Eigen_Dynamic.hpp"
#endif


#endif /* PROJECTS_LINALG_SRC_CL_MATRIX_HPP_ */
