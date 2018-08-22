/*
 * cl_XTK_Matrix_Teuchos.hpp
 *
 *  Created on: Aug 2, 2018
 *      Author: ktdoble
 */

#ifndef SRC_LINALG_CL_XTK_MATRIX_TEUCHOS_HPP_
#define SRC_LINALG_CL_XTK_MATRIX_TEUCHOS_HPP_

#ifdef XTK_USE_TEUCHOS
#include "Teuchos_SerialDenseMatrix.hpp"  // for SerialDenseMatrix^M
#include <climits>


namespace xtk
{

// NOTE: it appears as of 01/10/18 that SerialDenseMatrix doesn't like the ordinal type being size_t due to ambiguity in one of the functions it calls during construction
typedef Teuchos::SerialDenseMatrix<xtk::lint, xtk::real> Default_Matrix_Real;
typedef Teuchos::SerialDenseMatrix<xtk::lint, xtk::size_t> Default_Matrix_Integer;
typedef Teuchos::SerialDenseMatrix<xtk::lint, xtk::sint> Default_Matrix_Integer_Sint;

template<typename ScalarType>
Matrix_Type transpose(Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix)
{
    Teuchos::SerialDenseMatrix<xtk::lint,ScalarType>  newMatrix(aMatrix, Teuchos::ETransp::TRANS);
    return newMatrix;
}


template<typename Matrix_Type>
xtk::real determinant(Matrix_Type const & aMatrix)
{
    XTK_ASSERT(aMatrix.numRows() > 0, "To compute determinant, matrix must have nonzero number of rows");
    XTK_ASSERT(aMatrix.numCols() > 0, "To compute determinant, matrix must have nonzero number of cols");

    XTK_ASSERT(aMatrix.numRows()  == aMatrix.numCols(), "To compute determinant the matrix must be square");

//    XTK_ASSERT(aMatrix.numRows() <= 3, "To compute determinant, matrix must be 2x2 or 3x3");

    xtk::real detVal = 0;

    switch (aMatrix.numRows())
    {
        case 1:
            detVal = aMatrix(0,0);
            break;
        case 2:
        {
            detVal = aMatrix(0,0)*aMatrix(1,1)-aMatrix(0,1)*aMatrix(1,0);
            break;
        }
        default:
            {

                detVal = 0;
                xtk::real sign = 1.0;

                for (xtk::size_t colInd = 0; colInd < static_cast<xtk::size_t>(aMatrix.numCols()); ++colInd)
                {
                    Matrix_Type extractedMatrix = ExtractAllButASingleRowAndColumn(aMatrix, static_cast<xtk::size_t>(0), colInd);
                    detVal+=sign*determinant(extractedMatrix)*aMatrix(0,colInd);
                    sign*=-1.0;
                }

                break;
            }

    }
    return detVal;
}



template<typename ScalarType >
Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> operator-(Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix1,
                      Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix2)
{
    Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> newMatrix(aMatrix1);
    newMatrix-=aMatrix2;
    return newMatrix;
}


Teuchos::SerialDenseMatrix<xtk::lint,real> operator+(Teuchos::SerialDenseMatrix<xtk::lint,real> const & aMatrix1,
                                                     Teuchos::SerialDenseMatrix<xtk::lint,real> const & aMatrix2)
{
    Teuchos::SerialDenseMatrix<xtk::lint,real> newMatrix(aMatrix1);
    newMatrix+=aMatrix2;
    return newMatrix;
}


template<typename Type, typename Matrix_Type >
auto operator+(Matrix_Base<Type,Matrix_Type> & aMatrix1,
               Matrix_Base<Type,Matrix_Type> & aMatrix2)
->decltype(aMatrix1.matrix_data()+aMatrix2.matrix_data())
{
    return aMatrix1.matrix_data()+aMatrix2.matrix_data();
}

template<typename ScalarType >
Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> operator*(Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix1,
                      Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix2)
{
    Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> newMatrix(aMatrix1.numRows(),aMatrix2.numCols());
    newMatrix.multiply(Teuchos::ETransp::NO_TRANS,Teuchos::ETransp::NO_TRANS, 1.0, aMatrix1, aMatrix2, 0.0);
    return newMatrix;
}

template<typename ScalarType >
Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> operator*(ScalarType const & aScalar,
                      Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> const & aMatrix)
{
    Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> newMatrix(aMatrix);
    newMatrix*=aScalar;
    return newMatrix;
}

template<typename ScalarType>
Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> component_wise_abs(Teuchos::SerialDenseMatrix<xtk::lint,ScalarType> & aMatrix)
{
    size_t tNumRows = aMatrix.numRows();
    size_t tNumCols = aMatrix.numCols();

    for(size_t i = 0; i<tNumRows; i++)
    {
        for(size_t j = 0; j<tNumCols; j++)
        {
            aMatrix(i, j) = std::abs( aMatrix(i, j));
        }
    }

    return aMatrix;
}


template<typename Matrix_Type, typename Ordinal_Type>
Matrix_Type ExtractAllButASingleRowAndColumn(const Matrix_Type& aMatrix, Ordinal_Type aRowIndex, Ordinal_Type aColIndex)
{
    Matrix_Type extractedMatrix(aMatrix.numRows()-1,aMatrix.numCols()-1);

    xtk::size_t ouputRowInd = 0;

    for (xtk::size_t rowInd = 0; rowInd < static_cast<xtk::size_t>(aMatrix.numRows()); ++rowInd )
    {

        xtk::size_t ouputColInd = 0;
        if (rowInd != aRowIndex)
        {

            for (xtk::size_t colInd = 0; colInd < static_cast<xtk::size_t>(aMatrix.numCols()); ++colInd)
            {
                if (colInd != aColIndex)
                {
                    extractedMatrix(ouputRowInd, ouputColInd) = aMatrix(rowInd,colInd);
                    ++ouputColInd;
                }
            }
            ++ouputRowInd;
        }
    }

    return extractedMatrix;
}

template<typename Type, typename Matrix_Type>
class Matrix_Default : public xtk::Matrix_Base<Type, Matrix_Type>
{
public:
    Matrix_Default(size_t aNumRows,
                   size_t aNumColumns) :
            mMatrix(aNumRows, aNumColumns)
    {
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    Matrix_Default(Matrix_Type aBaseMatrix)
    {
        mMatrix = aBaseMatrix;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    void resize(size_t aNumRows, size_t aNumColumns)
    {
        mMatrix.reshape(aNumRows, aNumColumns);
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    void set_size(size_t aNumRows, size_t aNumColumns)
    {
        mMatrix.reshape(aNumRows, aNumColumns);
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    void fill(Type aFillValue)
    {
        mMatrix.putScalar(aFillValue);
    }

    // -----------------------------------------------------------------------------------------------------------------------------------------------
    //FIXME
    void set_row(size_t aRowIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aRow)
    {
        XTK_ASSERT(aRow.get_num_rows() == 1, "aRow needs to be a row matrix");
        XTK_ASSERT(aRowIndex < this->get_num_rows(), "Specified row index out of bounds");
        XTK_ASSERT(aRow.get_num_columns() == this->get_num_columns(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of columns)");

        size_t tROW_INDEX = 0;
//        mMatrix.row(aRowIndex) = dynamic_cast<const xtk::Matrix_Default<Type, Matrix_Type> &>(aRow).mMatrix.row(tROW_INDEX);

        for (xtk::size_t colIndex = 0; colIndex < this->get_num_columns(); ++colIndex)
        {
            mMatrix(aRowIndex,colIndex) = aRow(tROW_INDEX,colIndex);
        }

    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // FIXME
    void set_column(size_t aColumnIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aColumn)
    {

        XTK_ASSERT(aColumn.get_num_columns() == 1, "aColumn needs to be a column matrix");
        XTK_ASSERT(aColumnIndex < this->get_num_columns(), "Specified column index out of bounds");
        XTK_ASSERT(aColumn.get_num_rows() == this->get_num_rows(),
                   "Dimension mismatch (argument matrix and member matrix do not have same number of rows)");

        size_t tCOLUMN_INDEX = 0;
//        mMatrix.col(aColumnIndex) = dynamic_cast<const xtk::Matrix_Default<Type, Matrix_Type> &>(aColumn).mMatrix.col(tCOLUMN_INDEX);
        for (xtk::size_t rowIndex = 0; rowIndex < this->get_num_rows(); ++rowIndex)
        {
            mMatrix(rowIndex, aColumnIndex) = aColumn(rowIndex, tCOLUMN_INDEX);
        }

    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    //FIXME
    void get_row(size_t aRowIndex, xtk::Matrix_Base<Type, Matrix_Type> & aRow) const
    {
        if(aRow.get_num_rows() != 1)
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = aRow needs to be a row matrix\n\n";
        }
        if(aRowIndex >= this->get_num_rows())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Specified row index out of bounds\n\n";
        }
        if(aRow.get_num_columns() != this->get_num_columns())
        {
            std::cerr
                    << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                    << "MESSAGE = Dimension mismatch (argument matrix and member matrix do not have same number of columns)\n\n";
        }
//        xtk::Matrix_Default<Type, Matrix_Type> & tMatrix = dynamic_cast<xtk::Matrix_Default<Type, Matrix_Type> &>(aRow);
        const size_t tROW_INDEX = 0;
//        tMatrix.mMatrix.row(tROW_INDEX) = mMatrix.row(aRowIndex);
        for (xtk::size_t colIndex = 0; colIndex < this->get_num_columns(); ++colIndex)
        {
             aRow(tROW_INDEX,colIndex) = mMatrix(aRowIndex,colIndex);
        }
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // FIXME
    void get_column(size_t aColumnIndex, xtk::Matrix_Base<Type, Matrix_Type> & aColumn) const
    {
        if(aColumn.get_num_columns() != 1)
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = aColumn needs to be a column matrix\n\n";
        }
        if(aColumnIndex >= this->get_num_columns())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Specified column index out of bounds\n\n";
        }
        if(aColumn.get_num_rows() != this->get_num_rows())
        {
            std::cerr << "\n\nXTK ERROR IN FILE = " << __FILE__ << "FUNCTION = " << __FUNCTION__ << "LINE = " << __LINE__
                      << "MESSAGE = Dimension mismatch (argument matrix and member matrix do not have same number of rows)\n\n";
        }

//        xtk::Matrix_Default<Type, Matrix_Type> & tMatrix = dynamic_cast<xtk::Matrix_Default<Type, Matrix_Type> &>(aColumn);
        const size_t tCOLUMN_INDEX = 0;
//        tMatrix.mMatrix.col(tCOLUMN_INDEX) = mMatrix.col(aColumnIndex);

        Matrix_Type & aColumnMatrix = aColumn.matrix_data();
        XTK_ASSERT(this->get_num_rows() == static_cast<xtk::size_t>(mMatrix.stride()), "Ensure underlying matrix stride matches number of rows");
        std::copy( mMatrix[aColumnIndex],  mMatrix[aColumnIndex]+mMatrix.stride(), aColumnMatrix[tCOLUMN_INDEX]);
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    size_t get_num_rows() const
    {
        size_t tNumRows = mMatrix.numRows();
        return tNumRows;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    size_t get_num_columns() const
    {
        size_t tNumColumns = mMatrix.numCols();
        return tNumColumns;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    size_t get_num_matrix_elements() const
    {
        size_t tNumMatrixElements = mMatrix.numRows()*mMatrix.numCols();
        return tNumMatrixElements;
    }

    // -----------------------------------------------------------------------------------------------------------------------------------------------
    const Type* data() const
    {
        return mMatrix.values();
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    Matrix_Type & matrix_data()
    {
        return mMatrix;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    const Matrix_Type & matrix_data() const
        {
            return mMatrix;
        }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // FIXME...
    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(size_t aNumRows = 0,
                                              size_t aNumColumns = 0,
                                              Type aFillValue = std::numeric_limits<Type>::max()) const
    {
        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy;
        if(aNumRows == 0 && aNumColumns == 0)
        {
            size_t tNumRows = this->get_num_rows();
            size_t tNumColumns = this->get_num_columns();
            tMatrixCopy.reset(new Matrix_Default<Type, Matrix_Type>(tNumRows, tNumColumns));
        }

        else if(aNumColumns == 0)
        {
            size_t tNumColumns = this->get_num_columns();
            tMatrixCopy.reset(new Matrix_Default<Type, Matrix_Type>(aNumRows, tNumColumns));
        }

        else if(aNumRows == 0)
        {
            size_t tNumRows = this->get_num_rows();
            tMatrixCopy.reset(new Matrix_Default<Type, Matrix_Type>(tNumRows, aNumColumns));
        }

        else
        {
            tMatrixCopy.reset(new Matrix_Default<Type, Matrix_Type>(aNumRows, aNumColumns));
        }
        tMatrixCopy->fill(aFillValue);

        return tMatrixCopy;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // FIXME
    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(std::initializer_list<std::initializer_list<Type> > const & list) const
    {
        XTK_ASSERT(list.size() > 0, "The initializer list is empty.");

        size_t i = 0;
        size_t j = 0;
        size_t tNumRows = list.size();
        size_t tNumColumns = list.begin()->size();

        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = std::make_shared<Matrix_Default<Type, Matrix_Type>>(tNumRows, tNumColumns);

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

    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // FIXME
    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(Matrix_Type aBaseMatrix) const
    {
        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = std::make_shared<Matrix_Default<Type, Matrix_Type>>(aBaseMatrix);
        return tMatrixCopy;
    }


    // -----------------------------------------------------------------------------------------------------------------------------------------------
    std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> copy() const
    {
        size_t tNumRows = this->get_num_rows();
        size_t tNumColumns = this->get_num_columns();

        std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> tMatrixCopy = std::make_shared<Matrix_Default<Type,Matrix_Type>>(tNumRows, tNumColumns);

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
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    Type get_max_value() const
    {
        XTK_ASSERT(this->get_num_rows() > 0, "To get max value of matrix, matrix must have nonzero size");
        XTK_ASSERT(this->get_num_columns() > 0, "To get max value of matrix, matrix must have nonzero size");

        Type tMaxValue = mMatrix(0,0);
        for (xtk::size_t rowIndex = 0; rowIndex < this->get_num_rows(); ++rowIndex)
        {
            for (xtk::size_t colIndex = 0; colIndex < this->get_num_columns(); ++colIndex)
            {
//                fprintf(stdout,"matrix(%i,%j) = %2.2f\n",static_cast<int>(rowIndex), static_cast<int>(colIndex), mMatrix(rowIndex, colIndex));
                if (mMatrix(rowIndex, colIndex) > tMaxValue)
                {
                    tMaxValue = mMatrix(rowIndex, colIndex);
                }
            }
        }
        return tMaxValue;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    Type get_min_value() const
    {
        XTK_ASSERT(this->get_num_rows() > 0, "To get max value of matrix, matrix must have nonzero size");
        XTK_ASSERT(this->get_num_columns() > 0, "To get max value of matrix, matrix must have nonzero size");

        Type tMinValue = mMatrix(0,0);
        for (xtk::size_t rowIndex = 0; rowIndex < this->get_num_rows(); ++rowIndex)
        {
            for (xtk::size_t colIndex = 0; colIndex < this->get_num_columns(); ++colIndex)
            {
                if (mMatrix(rowIndex, colIndex) < tMinValue)
                {
                    tMinValue = mMatrix(rowIndex, colIndex);
                }
            }
        }
        return tMinValue;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    Type & operator()(size_t aRowIndex, size_t aColumnIndex)
    {
        XTK_ASSERT(aRowIndex < this->get_num_rows(), "Requested row is out of bounds");
        XTK_ASSERT(aColumnIndex < this->get_num_columns(), "Requested column is out of bounds");
        Type & tValueReference = mMatrix(aRowIndex, aColumnIndex);
        return tValueReference;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    const Type & operator()(size_t aRowIndex, size_t aColumnIndex) const
    {
        XTK_ASSERT(aRowIndex < this->get_num_rows(), "Requested row is out of bounds");
        XTK_ASSERT(aColumnIndex < this->get_num_columns(), "Requested column is out of bounds");
        const Type & tValueReference = mMatrix(aRowIndex, aColumnIndex);
        return tValueReference;
    }
    // -----------------------------------------------------------------------------------------------------------------------------------------------
    // Private Member variables
    Matrix_Type mMatrix;
private:

    // Private Functions
private:
    Matrix_Default(const xtk::Matrix_Default<Type,Matrix_Type> &);  // copy constructor
    xtk::Matrix_Default<Type, Matrix_Type> & operator=(const xtk::Matrix_Default<Type,Matrix_Type> &); // copy operator
};

#endif

#endif /* SRC_LINALG_CL_XTK_MATRIX_TEUCHOS_HPP_ */
