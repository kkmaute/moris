
/*
 * cl_Linalg_Matrix.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */

#ifndef INCLUDE_CL_LINALG_MATRIX_HPP_
#define INCLUDE_CL_LINALG_MATRIX_HPP_

#include <memory> // for unique_ptr
#include <initializer_list>

namespace moris
{
template<typename Type, typename Matrix_Type>
class Matrix_Base
{
public:


    virtual ~Matrix_Base()
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
    virtual
    void
    resize( moris::size_t const & aNumRows,
            moris::size_t const & aNumCols ) = 0;

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
    virtual
    void set_size( moris::size_t aNumRows,
                   moris::size_t aNumColumns ) = 0;

    // -------------------------------------------------------------------------

    /**
     * Initialization operator.
     *
     * @param[in] aVal A given value to set all elements to that.
     *
     * Example:
     * @include cl_Mat/Mat_fill.inc
     */
    virtual
    void
    fill( Type const & aFillValue ) = 0;

    // -------------------------------------------------------------------------

    virtual void
    set_row(size_t aRowIndex,
            const moris::Matrix_Base<Type, Matrix_Type> & aRow) = 0;

    // -------------------------------------------------------------------------

    virtual void set_column(size_t aColumnIndex, const moris::Matrix_Base<Type, Matrix_Type> & aColumn) = 0;

    // -------------------------------------------------------------------------

    /**
     * Get row vector, similar to Matlab A(aK,:).
     *
     * @param[in] aK row index.
     *
     * @return row vector.
     *
     *
     *  Note: this is not a moris mat return
     */
    virtual auto row(moris::size_t const & aRowIndex) const = 0;

    // -------------------------------------------------------------------------

    virtual void get_column(size_t aColumnIndex, moris::Matrix_Base<Type, Matrix_Type> & aColumn) const = 0;

    // -------------------------------------------------------------------------
    /**
     * Get the number of columns in a data set, similar to Matlab cols().
     *
     * @return Number of columns.
     */
    virtual
    moris::size_t
    n_cols() const = 0;


    // -------------------------------------------------------------------------

     /**
     * Returns the number of elements in the %matrix.
     *
     * @return Number of elements in the %matrix.
     *
     */
    virtual
    size_t
    numel() const = 0;

    // -------------------------------------------------------------------------

    virtual
    const Type*
    data() const = 0;

    // -------------------------------------------------------------------------

    virtual
    Matrix_Type &
    matrix_data() = 0;

    // -------------------------------------------------------------------------

    virtual
    std::shared_ptr<moris::Matrix_Base<Type, Matrix_Type>>
    create(size_t aNumRows = 0, size_t aNumColumns = 0, Type aFillValue = 0) const = 0;

    // -------------------------------------------------------------------------

    virtual
    std::shared_ptr<moris::Matrix_Base<Type, Matrix_Type>>
    create(std::initializer_list< std::initializer_list< Type > > const & list) const = 0;

    // -------------------------------------------------------------------------

    virtual
    std::shared_ptr<moris::Matrix_Base<Type, Matrix_Type>>
    create(Matrix_Type aBaseMatrix) const = 0;

    // -------------------------------------------------------------------------

     virtual
     std::shared_ptr<moris::Matrix_Base<Type, Matrix_Type>>
     copy() const  = 0;

    // -------------------------------------------------------------------------

    /**
     * Member function.
     *
     * @return The extremum value of an object.
     *
     * Example:
     * @include cl_Mat/Mat_max.inc
     */
    virtual
    Type
    max() const = 0;

    // -------------------------------------------------------------------------

    /**
     * @param[in] aRowIndex Row index of max value of an object.
     * @param[in] aColIndex Column index of max value of an object.
     *
     * @return The extremum value of an object and store the location
     * of the extremum value in the provided variable(s)
     *
     */

    virtual
    Type
    max( moris::uint & aRowIndex,
         moris::uint & aColIndex )  const = 0;

    // -------------------------------------------------------------------------

    /**
     * Member function.
     *
     * @return The extremum value of an object.
     *
     * Example:
     * @include cl_Mat/Mat_min.inc
     */
    virtual
    Type
    min() const = 0;

    /**
     * Member function.
     * @param[in] aRowIndex Row index of min value of an object.
     * @param[in] aColIndex Column index of min value of an object.
     *
     * @return The extremum value of an object and store the location
     * of the extremum value in the provided variable(s)
     *
     */
    virtual
    Type
    min( moris::uint & aRowIndex,
         moris::uint & aColIndex ) const = 0;

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] i_index Row index for which data should be accessed.
     * @param[in] j_index Column index for which data should be accessed.
     */

    virtual
    Type &
    operator()( moris::size_t const & i_index,
                moris::size_t const & j_index ) = 0;

    // -------------------------------------------------------------------------

    /**
     * @brief Overloaded moris::Matrix_Base::operator()
     *
     * @param[in] i_index Row index for which data should be accessed.
     * @param[in] j_index Column index for which data should be accessed.
     */

    virtual
    const Type &
    operator()( moris::size_t const & i_index,
                moris::size_t const & j_index ) const = 0;

};
}


#endif /* INCLUDE_CL_LINALG_MATRIX_HPP_ */
