/*
 * cl_XTK_Matrix_Eigen.hpp
 *
 *  Created on: Jun 12, 2017
 *      Author: ktdoble
 */

#ifndef SRC_LINALG_CL_XTK_MATRIX_EIGEN_HPP_
#define SRC_LINALG_CL_XTK_MATRIX_EIGEN_HPP_
// Standard headers
#include <memory>   // Shared ptrs
#include <initializer_list>
#include <type_traits>

// Eigen headers
#include <Eigen/Dense>

#include "assert/fn_xtk_assert.hpp"
#include "linalg/cl_XTK_Matrix_Base.hpp"
#include "cl_Matrix.hpp"

namespace xtk
{
typedef moris::DDRMat   moris::DDRMat;
typedef moris::DDSTMat  moris::DDSTMat;
typedef moris::DDSMat   moris::DDSTMat_Sint;


/*
 * EIGEN SPECIFIC FREE FUNCTIONS
 */
template<typename Eigen_Exp>
auto transpose(Eigen_Exp  const & aMatrix)
->decltype(aMatrix.transpose())
{
    return aMatrix.transpose();
}

template<typename Eigen_Exp>
auto determinant(Eigen_Exp const & aMatrix)
->decltype(aMatrix.determinant())
{
    return aMatrix.determinant();
}


template<typename Eigen_Exp>
auto component_wise_abs(Eigen_Exp & aExpTemplate)
->decltype(aExpTemplate.cwiseAbs())
{
    return aExpTemplate.cwiseAbs();
}


template<typename Matrix_Type>
class Matrix_Base_Eigen : public xtk::Matrix_Base<Type, Matrix_Type>
{
public:
    Matrix_Base_Eigen(size_t const & aNumRows,
                      size_t const & aNumColumns) :
            mMatrix(aNumRows, aNumColumns)
    {
    }

    Matrix_Base_Eigen(Matrix_Type aBaseMatrix)
    {
        mMatrix = aBaseMatrix;
    }

    Matrix_Base_Eigen(std::initializer_list<std::initializer_list<Type> > const & aInitList)
    {
        size_t i = 0;
        size_t j = 0;
        size_t tNumRows = aInitList.size();
        size_t tNumColumns = aInitList.begin()->size();

        mMatrix = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>(tNumRows, tNumColumns);

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


    void resize(size_t aNumRows, size_t aNumColumns)
    {
        mMatrix.conservativeResize(aNumRows, aNumColumns);
    }

    void set_size(size_t aNumRows, size_t aNumColumns)
    {
        mMatrix.resize(aNumRows, aNumColumns);
    }

    void fill(Type aFillValue)
    {
        mMatrix.fill(aFillValue);
    }

    void set_row(size_t aRowIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aRow)
    {
        XTK_ASSERT(aRow.n_rows() == 1, "aRow needs to be a row matrix");
        XTK_ASSERT(aRowIndex < this->n_rows(), "Specified row index out of bounds");
        XTK_ASSERT(aRow.n_cols() == this->n_cols(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

        size_t tROW_INDEX = 0;
        mMatrix.row(aRowIndex) = dynamic_cast<const xtk::Matrix_Base_Eigen<Type, Matrix_Type> &>(aRow).mMatrix.row(tROW_INDEX);
    }

    void set_column(size_t aColumnIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aColumn)
    {

        XTK_ASSERT(aColumn.n_cols() == 1, "aColumn needs to be a column matrix");
        XTK_ASSERT(aColumnIndex < this->n_cols(), "Specified column index out of bounds");
        XTK_ASSERT(aColumn.n_rows() == this->n_rows(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of rows)");

        size_t tCOLUMN_INDEX = 0;
        mMatrix.col(aColumnIndex) = dynamic_cast<const xtk::Matrix_Base_Eigen<Type, Matrix_Type> &>(aColumn).mMatrix.col(tCOLUMN_INDEX);
    }

    void get_row(size_t aRowIndex, xtk::Matrix_Base<Type, Matrix_Type> & aRow) const
    {
        if(aRow.n_rows() != 1)
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = aRow needs to be a row matrix\n\n";
        }
        if(aRowIndex >= this->n_rows())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Specified row index out of bounds\n\n";
        }
        if(aRow.n_cols() != this->n_cols())
        {
            std::cerr
                    << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                    << "MESSAGE = Dimension mismatch (argument matrix and member matrix do not have same number of columns)\n\n";
        }
        xtk::Matrix_Base_Eigen<Type, Matrix_Type> & tMatrix = dynamic_cast<xtk::Matrix_Base_Eigen<Type, Matrix_Type> &>(aRow);
        const size_t tROW_INDEX = 0;
        tMatrix.mMatrix.row(tROW_INDEX) = mMatrix.row(aRowIndex);
    }

    void get_column(size_t aColumnIndex, xtk::Matrix_Base<Type, Matrix_Type> & aColumn) const
    {
        if(aColumn.n_cols() != 1)
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = aColumn needs to be a column matrix\n\n";
        }
        if(aColumnIndex >= this->n_cols())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Specified column index out of bounds\n\n";
        }
        if(aColumn.n_rows() != this->n_rows())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Dimension mismatch (argument matrix and member matrix do not have same number of rows)\n\n";
        }

        xtk::Matrix_Base_Eigen<Type, Matrix_Type> & tMatrix = dynamic_cast<xtk::Matrix_Base_Eigen<Type, Matrix_Type> &>(aColumn);
        const size_t tCOLUMN_INDEX = 0;
        tMatrix.mMatrix.col(tCOLUMN_INDEX) = mMatrix.col(aColumnIndex);
    }

    size_t n_rows() const
    {
        size_t tNumRows = mMatrix.rows();
        return tNumRows;
    }

    size_t n_cols() const
    {
        size_t tNumColumns = mMatrix.cols();
        return tNumColumns;
    }

    size_t get_num_matrix_elements() const
    {
        size_t tNumMatrixElements = mMatrix.size();
        return tNumMatrixElements;
    }

    const Type* data() const
    {
        return mMatrix.data();
    }

    Matrix_Type & matrix_data()
    {
        return mMatrix;
    }

    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(size_t aNumRows = 0,
                                              size_t aNumColumns = 0,
                                              Type aFillValue = std::numeric_limits<Type>::max()) const
    {
        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy;
        if(aNumRows == 0 && aNumColumns == 0)
        {
            size_t tNumRows = this->n_rows();
            size_t tNumColumns = this->n_cols();
            tMatrixCopy.reset(new Matrix_Base_Eigen<Type, Matrix_Type>(tNumRows, tNumColumns));
        }

        else if(aNumColumns == 0)
        {
            size_t tNumColumns = this->n_cols();
            tMatrixCopy.reset(new Matrix_Base_Eigen<Type, Matrix_Type>(aNumRows, tNumColumns));
        }

        else if(aNumRows == 0)
        {
            size_t tNumRows = this->n_rows();
            tMatrixCopy.reset(new Matrix_Base_Eigen<Type, Matrix_Type>(tNumRows, aNumColumns));
        }

        else
        {
            tMatrixCopy.reset(new Matrix_Base_Eigen<Type, Matrix_Type>(aNumRows, aNumColumns));
        }
        tMatrixCopy->fill(aFillValue);

        return tMatrixCopy;
    }

    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(std::initializer_list<std::initializer_list<Type> > const & list) const
    {
        XTK_ASSERT(list.size() > 0, "The initializer list is empty.");

        size_t i = 0;
        size_t j = 0;
        size_t tNumRows = list.size();
        size_t tNumColumns = list.begin()->size();

        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = std::make_shared<Matrix_Base_Eigen<Type, Matrix_Type>>(tNumRows, tNumColumns);

        for(const auto tRow : list) // loop over number of rows
        {
            XTK_ASSERT(tRow.size() == list.begin()->size(),
                       "The number of elements in one of the rows does not equal the number of columns.");

            for(const auto tCol : tRow) // loop over every value in the row
            {
                (*tMatrixCopy)(i, j) = tCol;
                ++j;
            }
            j = 0;
            ++i;
        }

        return tMatrixCopy;
    }

    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(Matrix_Type aBaseMatrix) const
    {
        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = std::make_shared<Matrix_Base_Eigen<Type, Matrix_Type>>(aBaseMatrix);
        return tMatrixCopy;
    }

    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> copy() const
    {
        size_t tNumRows = this->n_rows();
        size_t tNumColumns = this->n_cols();

        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = this->create(tNumRows, tNumColumns);

        // Fill
        for( size_t r = 0; r<tNumRows; r++ )
        {
            for( size_t c = 0; c<tNumColumns; c++ )
            {
                (*tMatrixCopy)(r,c) = mMatrix(r,c);
            }
        }
        return tMatrixCopy;
    }

    Type get_max_value() const
    {
        Type tMaxValue = mMatrix.maxCoeff();
        return tMaxValue;
    }

    Type get_min_value() const
    {
        Type tMinValue = mMatrix.minCoeff();
        return tMinValue;
    }

    Type & operator()(size_t aRowIndex, size_t aColumnIndex)
    {
        XTK_ASSERT(aRowIndex < this->n_rows(), "Requested row is out of bounds");
        XTK_ASSERT(aColumnIndex < this->n_cols(), "Requested column is out of bounds");
        Type & tValueReference = mMatrix(aRowIndex, aColumnIndex);
        return tValueReference;
    }

    const Type & operator()(size_t aRowIndex, size_t aColumnIndex) const
    {
        XTK_ASSERT(aRowIndex < this->n_rows(), "Requested row is out of bounds");
        XTK_ASSERT(aColumnIndex < this->n_cols(), "Requested column is out of bounds");
        const Type & tValueReference = mMatrix(aRowIndex, aColumnIndex);
        return tValueReference;
    }

    // Private Member variables
private:
    Matrix_Type mMatrix;

    // Private Functions
private:
    Matrix_Base_Eigen(const xtk::Matrix_Base_Eigen<Type,Matrix_Type> &);  // copy constructor
    xtk::Matrix_Base_Eigen<Type, Matrix_Type> & operator=(const xtk::Matrix_Base_Eigen<Type,Matrix_Type> &); // copy operator
};
}

#endif /* SRC_LINALG_CL_XTK_MATRIX_EIGEN_HPP_ */
