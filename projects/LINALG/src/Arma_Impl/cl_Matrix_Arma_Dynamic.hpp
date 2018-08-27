/*
 * cl_Matrix_Arma_3x3.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_ARMA_IMPL_CL_MATRIX_ARMA_DYNAMIC_HPP_
#define PROJECTS_LINALG_SRC_ARMA_IMPL_CL_MATRIX_ARMA_DYNAMIC_HPP_
#include <armadillo>

#include "typedefs.hpp"

#include "cl_Matrix.hpp"

namespace moris
{

template<typename Type>
class Mat_New<Type, arma::Mat<Type>>
{
private:
    arma::Mat<Type> mMatrix;
public:
    Mat_New(){}

    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols):
                mMatrix(aNumRows,aNumCols)
    {

    }

    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols,
            Type   const & aFillVal):
            mMatrix( aNumRows, aNumCols, aFillVal )
    {

    }

    // template constructor
    template< typename A >
    Mat_New(A const & X ):
            mMatrix(X)
    {

    }

    Mat_New(std::initializer_list<std::initializer_list<Type> > const & aInitList)
    {
        size_t i = 0;
        size_t j = 0;
        size_t tNumRows = aInitList.size();
        size_t tNumColumns = aInitList.begin()->size();

        mMatrix = arma::Mat<Type>(tNumRows, tNumColumns);

        for(const auto tRow : aInitList) // loop over number of rows
        {
            XTK_ASSERT(tRow.size() == aInitList.begin()->size(),
                       "The number of elements in one of the rows does not equal the number of columns.");

            for(const auto tCol : tRow) // loop over every value in the row
            {
                mMatrix(i, j) = tCol;
                ++j;
            }
            j = 0;
            ++i;
        }
    }

    void
    resize(const size_t & aNumRows,
           const size_t & aNumCols)
    {
        mMatrix.conservativeResize(aNumRows, aNumCols);
    }

    void
    fill(const Type & aFillValue)
    {
        mMatrix.fill(aFillValue);
    }

    /**
     * Get the number of columns in a data set, similar to Matlab cols().
     *
     * @return Number of columns.
     */
    size_t
    n_cols()
    {
        return mMatrix.n_cols;
    }

    /**
     * Get the number of rows in a data set, similar to Matlab rows().
     *
     * @return Number of rows.
     */
    size_t
    n_rows()
    {
        return mMatrix.n_rows;
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
        return mMatrix.n_elem;
    }

    void set_row(size_t aRowIndex, const Mat_New<Type, arma::Mat<Type>> & aRow)
    {
        MORIS_ASSERT(aRow.get_num_rows() == 1, "aRow needs to be a row matrix");
        MORIS_ASSERT(aRowIndex < this->get_num_rows(), "Specified row index out of bounds");
        MORIS_ASSERT(aRow.get_num_columns() == this->get_num_columns(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

        size_t tROW_INDEX = 0;
        mMatrix.row(aRowIndex) = aRow.matrix_data().row(tROW_INDEX);
    }

    void set_column(size_t aColumnIndex, const Mat_New<Type, arma::Mat<Type>> & aColumn)
    {

        MORIS_ASSERT(aColumn.get_num_columns() == 1, "aColumn needs to be a column matrix");
        MORIS_ASSERT(aColumnIndex < this->get_num_columns(), "Specified column index out of bounds");
        MORIS_ASSERT(aColumn.get_num_rows() == this->get_num_rows(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of rows)");

        size_t tCOLUMN_INDEX = 0;
        mMatrix.col(aColumnIndex) = aColumn.matrix_data().col(tCOLUMN_INDEX);
    }


    void
    get_column(size_t aColumnIndex,
               Mat_New<Type, arma::Mat<Type>> & aColumn) const
    {
        MORIS_ASSERT(aColumn.n_cols() == 1,"aColumn needs to be a column matrix");
        MORIS_ASSERT(aColumnIndex < this->n_cols(),"Specified column index out of bounds");
        MORIS_ASSERT(aColumn.n_rows() == this->n_rows(),"Dimension mismatch (argument matrix and member matrix do not have same number of rows)");
        const size_t tCOLUMN_INDEX = 0;
        aColumn.matrix_data().col(tCOLUMN_INDEX) = mMatrix.col(aColumnIndex);
    }

    void get_row(size_t aRowIndex, Mat_New<Type, arma::Mat<Type>> & aRow) const
    {
        MORIS_ASSERT(aRow.n_rows() == 1,"aRow needs to be a row matrix");
        MORIS_ASSERT(aRowIndex < this->n_rows(),"Specified row index out of bounds");
        MORIS_ASSERT(aRow.n_cols() == this->n_cols(),"Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

        const size_t tROW_INDEX = 0;
        aRow.mMatrix.row(tROW_INDEX) = mMatrix.row(aRowIndex);
    }

    Type*
    data() const
    {
        return mMatrix.data();
    }

    inline
    arma::Mat<Type> &
    matrix_data()
    {
        return mMatrix;
    }

    Type
    max() const
    {
        MORIS_ASSERT(false,"Entered non-specialized base class of Matrix, Has your matrix_type template been implemented and the correct header included?");
        return 0;
    }

    inline
    Type &
    operator()( size_t const & i_index,
                size_t const & j_index )
    {
        return mMatrix(i_index,j_index);
    }


};
}





#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_CL_MATRIX_ARMA_DYNAMIC_HPP_ */
