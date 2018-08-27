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
template<typename Type, typename Matrix_Type>
class Mat_New
{
    // These member variables are here for error throwing
    Matrix_Type mMatrix;
    Type*       mDummy;

public:
    Mat_New()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // Constructor with no fill value
    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // Constructor with fill value
    Mat_New(size_t const & aNumRows,
        size_t const & aNumCols,
        Type   const & aFillVal)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    /* Constructor for std::initializer_list<std::initializer_list>,
     * allows for {{1,2,3},{4,5,6},{7,8,9}} to be passed as input
     * and receive,
     * [ 1 2 3 ]
     * [ 4 5 6 ]
     * [ 7 8 9 ]
     *
     * as output
    */
    Mat_New(std::initializer_list<std::initializer_list<Type> > const & aInitList)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");

    }

    // template constructor
    template< typename A >
    Mat_New(A const & X ):
                mMatrix(X)
     {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");

     }

    // Copy operations
    Mat_New<Type,Matrix_Type>
    copy()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
    resize(const size_t & aNumRows,
           const size_t & aNumCols)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
    void set_size( size_t aNumRows,
                   size_t aNumColumns )
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
    fill(const Type & aFillValue)
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // -------------------------------------------------------------------------

    void
    set_row( size_t aRowIndex,
             const Mat_New<Type, Matrix_Type> & aRow )
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // -------------------------------------------------------------------------

    void
    set_column( size_t                             aColumnIndex,
                const Mat_New<Type, Matrix_Type> & aColumn )
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // -------------------------------------------------------------------------

    void
    get_column(size_t aColumnIndex,
               Mat_New<Type, Matrix_Type> & aColumn) const
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
    }

    // -------------------------------------------------------------------------
    /**
     * Get the number of columns in a data set, similar to Matlab cols().
     *
     * @return Number of columns.
     */
    size_t
    n_cols()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return 0;
    }

    /**
     * Get the number of rows in a data set, similar to Matlab rows().
     *
     * @return Number of rows.
     */
    size_t
    n_rows()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return 0;
    }

    // -------------------------------------------------------------------------

    Type*
    data() const
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return mDummy;
    }

    // -------------------------------------------------------------------------

    Matrix_Type &
    matrix_data()
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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

    Type
    max() const
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return 0;
    }

    // -------------------------------------------------------------------------

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
    // -------------------------------------------------------------------------

    /**
     * Member function.
     *
     * @return The extremum value of an object.
     *
     * Example:
     * @include cl_Mat/Mat_min.inc
     */
    Type
    min() const
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return 0;
    }
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
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
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
