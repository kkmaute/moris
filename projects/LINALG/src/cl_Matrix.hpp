/*
 * cl_Matrix.hpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_CL_MATRIX_HPP_
#define PROJECTS_LINALG_SRC_CL_MATRIX_HPP_

#include "cl_Matrix_Base.hpp"
#include "cl_Dense_Backend_Matrix_Factory.hpp"

#include <memory>
#include "typedefs.hpp"

namespace moris
{

template<typename Type, typename Matrix_Type>
class Mat_New
{
public:
    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols)
    {
        A<Type, Matrix_Type> tA;
//        Dense_Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
        mBaseMat = tA.create_matrix_base(aNumRows,aNumCols);
    }

    std::shared_ptr<Matrix_Base<Type,Matrix_Type>> mBaseMat;

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
    resize(const moris::size_t & aNumRows,
           const moris::size_t & aNumCols)
    {
        mBaseMat->resize(aNumRows,aNumCols);
    }

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
//    void set_size( moris::size_t aNumRows,
//                   moris::size_t aNumColumns )
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
//    moris::size_t
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
//    moris::size_t
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
//    Matrix_Type &
//    matrix_data()
//    {
//        return mBaseMat->matrix_data();
//    }
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
    operator()(const moris::size_t & aRowIndex,
               const moris::size_t & aColIndex)
    {
        return (*mBaseMat)(aRowIndex,aColIndex);
    }

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] aRowIndex Row index for which data should be accessed.
     * @param[in] aColIndex Column index for which data should be accessed.
     */

    const Type &
    operator()(const moris::size_t & aRowIndex,
               const moris::size_t & aColIndex) const
    {
        return (*mBaseMat)(aRowIndex,aColIndex);

    }

};
}
#endif /* PROJECTS_LINALG_SRC_CL_MATRIX_HPP_ */
