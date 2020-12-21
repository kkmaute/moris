/*
 * cl_XTK_Matrix.hpp
 *
 *  Created on: Aug 7, 2017
 *      Author: ktdoble
 */

#ifndef SRC_LINALG_CL_XTK_MATRIX_HPP_
#define SRC_LINALG_CL_XTK_MATRIX_HPP_

#include "cl_Matrix.hpp"

//#include <memory> // for unique_ptr
//#include "core/cl_XTK_Parameters.hpp"
//#include "cl_Matrix.hpp"
//#include "linalg/cl_XTK_Backend_Matrix_Factory.hpp"
//#include "cl_XTK_Enums.hpp"

//namespace xtk
//{
//template < typename Type,
//           typename Matrix_Type >
//class Mat
//{
//public:
//
//    Mat(){};
//
//    Mat(size_t const & aNumRows,
//        size_t const & aNumCols)
//    {
//        Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
//        mMatrixPointer = tBEFactory.create_matrix_base(aNumRows,aNumCols);
//    }
//
//    Mat(size_t const & aNumRows,
//        size_t const & aNumCols,
//        Type   const & aFilVal)
//    {
//        Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
//        mMatrixPointer = tBEFactory.create_matrix_base(aNumRows,aNumCols);
//        fill(aFilVal);
//    }
//
//    Mat(std::initializer_list<std::initializer_list<Type> > const & aInitList)
//    {
//        Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
//        mMatrixPointer = tBEFactory.create_matrix_base(aInitList);
//    }
//
//    // Explicit to prevent copies unless desired to
//    explicit
//    Mat(Matrix_Base<Type,Matrix_Type> const & aMatBase)
//    {
//        size_t tNumRow = aMatBase.n_rows();
//        size_t tNumCol = aMatBase.n_cols();
//        Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
//        mMatrixPointer = tBEFactory.create_matrix_base(tNumRow,tNumCol);
//        for(size_t i = 0; i<tNumRow; i++)
//        {
//            for(size_t j = 0; j<tNumCol; j++)
//            {
//                (*mMatrixPointer)(i,j) = aMatBase(i,j);
//            }
//        }
//    }
//
//    explicit
//    Mat(Matrix_Type const & aMatBase)
//    {
//        Backend_Matrix_Factory<Type,Matrix_Type> tBEFactory;
//        mMatrixPointer = tBEFactory.create_matrix_base(aMatBase);
//    }
//
//
//    ~Mat()
//    {
//    }
//
//    //Required Container like properties
//
//    /*
//     * Resizes the matrix and preserves the data
//     * @param[in] aNumRows - Number of rows
//     * @param[in] aNumColumns - Number of columns
//     */
//    void
//    resize(size_t aNumRows, size_t aNumColumns)
//    {
//        mMatrixPointer->resize(aNumRows,aNumColumns);
//    }
//
//    /*
//     * Resizes the matrix and does not necessarily preserves the data
//     * @param[in] aNumRows - Number of rows
//     * @param[in] aNumColumns - Number of columns
//     */
//    void set_size(size_t aNumRows, size_t aNumColumns)
//    {
//        mMatrixPointer->set_size(aNumRows,aNumColumns);
//    }
//
//    /*
//     * Fill the matrix with a value
//     * @param[in] aFillValue - value to fill with
//     */
//    void fill(Type aFillValue)
//    {
//        mMatrixPointer->fill(aFillValue);
//    }
//
//    /*
//     * Replaces row in the matrix
//     * @param[in] aRowIndex - Row to set
//     * @param[in] aRow - Data to set in row
//     */
//    void set_row(size_t aRowIndex,
//                 Mat<Type, Matrix_Type> const & aRow)
//    {
//        mMatrixPointer->set_row(aRowIndex,aRow.matrix_base());
//    }
//
//    /*
//     * Replaces column in matrix
//     * @param[in] aColumnIndex - Column to set
//     * @param[in] aColumn - Data to set in row
//     */
//    void set_column(size_t aColumnIndex,
//                    Mat<Type, Matrix_Type> const & aColumn)
//    {
//        mMatrixPointer->set_column(aColumnIndex,aColumn.matrix_base());
//    }
//
//    /*
//     * Returns a copy of row
//     * @param[in] aRowIndex - Row index of matrix to get
//     * @param[out] aRow - Row to Copy into.
//     */
//    Mat<Type,Matrix_Type>
//    get_row(size_t aRowIndex) const
//    {
//        Mat<Type,Matrix_Type> tRow(1,this->n_cols());
//
//        mMatrixPointer->get_row(aRowIndex, tRow.matrix_base());
//
//        return tRow;
//    }
//
//    /*
//     * Returns a copy of column
//     * @param[in] aColumnIndex - Column index of matrix to get
//     * @param[out] aColumn -  Column to copy into
//     */
//    Mat<Type,Matrix_Type>
//    get_column(size_t aColumnIndex) const
//    {
//        Mat<Type,Matrix_Type> tCol(this->n_rows(),1);
//
//        mMatrixPointer->get_column(aColumnIndex, tCol.matrix_base());
//
//        return tCol;
//    }
//
//    /*
//     * Returns the number of rows in matrix
//     */
//    size_t n_rows() const
//    {
//        return mMatrixPointer->n_rows();
//    }
//
//    /*
//     * Returns the number of columns in matrix
//     */
//    size_t n_cols() const
//    {
//        return mMatrixPointer->n_cols();
//    }
//
//    /*
//     * Returns the number of elements in the matrix
//     */
//    size_t get_num_elements() const
//    {
//        return mMatrixPointer->get_num_matrix_elements();
//    }
//
//    /*
//     * Returns the maximum value in the matrix
//     */
//    Type get_max_value() const
//    {
//        return mMatrixPointer->get_max_value();
//    }
//
//    /*
//     * Returns the minimum value in the matrix
//     */
//    Type get_min_value() const
//    {
//        return mMatrixPointer->get_min_value();
//    }
//
//    /*
//     * Access the data at a location non const
//     * @param[in] aRowIndex - Row index,
//     * @param[in] aColumnsIndex - Column index
//     */
//    Mat<Type,Matrix_Type> operator = (Matrix_Type const & aBackendMat)
//    {
//        return Mat<Type,Matrix_Type>(aBackendMat);
//    }
//
//    /*
//     * Access the data at a location non const
//     * @param[in] aRowIndex - Row index,
//     * @param[in] aColumnsIndex - Column index
//     */
//    Type & operator()(size_t aRowIndex, size_t aColumnIndex)
//    {
//        return (*mMatrixPointer)(aRowIndex,aColumnIndex);
//    }
//
//    /*
//     * Access the data at a location const
//     * @param[in] aRowIndex - Row index,
//     * @param[in] aColumnsIndex - Column index
//     */
//    const Type & operator()(size_t aRowIndex, size_t aColumnIndex) const
//    {
//        return (*mMatrixPointer)(aRowIndex,aColumnIndex);
//    }
//
//
//    /*
//     * Makes a copy of matrix
//     */
//    Mat<Type,Matrix_Type>
//    copy() const
//    {
//        Mat<Type,Matrix_Type> tMatCopy(this->n_rows(),
//                                       this->n_cols());
//
//        for( size_t i = 0 ; i <this->n_rows(); i++)
//        {
//            for( size_t j = 0 ; j <this->n_cols(); j++)
//            {
//                tMatCopy(i,j) = (*mMatrixPointer)(i,j);
//            }
//        }
//
//        return tMatCopy;
//    }
//
//
//    /*
//     * Returns the underlying array of data
//     */
//    const Type* data() const
//    {
//        return mMatrixPointer->data();
//    }
//
//    /*
//     * Returns the Matrix Base
//     */
//    Matrix_Base<Type,Matrix_Type> &
//    matrix_base()
//    {
//        return *mMatrixPointer;
//    }
//
//    Matrix_Base<Type,Matrix_Type> const &
//    matrix_base() const
//    {
//        return *mMatrixPointer;
//    }
//
//    /*
//     * Returns the underlying matrix (eigen/teuchos)
//     */
//    Matrix_Type &  matrix_data()
//    {
//        return mMatrixPointer->matrix_data();
//    }
//
//    std::shared_ptr<Matrix_Base<Type,Matrix_Type>> mMatrixPointer;
//
//};
//
//template<typename Type, typename Matrix_Type >
//auto operator-(Mat<Type,Matrix_Type> const & aMatrix1,
//               Mat<Type,Matrix_Type> const & aMatrix2)
//->decltype(aMatrix1.matrix_data()-aMatrix2.matrix_data())
//{
//    return aMatrix1.matrix_data()-aMatrix2.matrix_data();
//}
//
//template<typename Type, typename Matrix_Type >
//auto operator+(Mat<Type,Matrix_Type> & aMatrix1,
//               Mat<Type,Matrix_Type> & aMatrix2)
//->decltype(aMatrix1.matrix_data()+aMatrix2.matrix_data())
//{
//    return aMatrix1.matrix_data()+aMatrix2.matrix_data();
//}
//
//
//template<typename Type, typename Matrix_Type >
//auto operator*(Mat<Type,Matrix_Type> & aMatrix1,
//               Mat<Type,Matrix_Type> & aMatrix2)
//->decltype(aMatrix1.matrix_data()+aMatrix2.matrix_data())
//{
//    return aMatrix1.matrix_data()+aMatrix2.matrix_data();
//}
//
//template<typename Type, typename Matrix_Type, typename Exp_Temp >
//auto operator*(Exp_Temp & aMatrix1,
//               Mat<Type,Matrix_Type> & aMatrix2)
//->decltype(aMatrix1+aMatrix2.matrix_data())
//{
//    return aMatrix1+aMatrix2.matrix_data();
//}
//
//template<typename Type, typename Matrix_Type, typename Exp_Temp >
//auto operator*(Mat<Type,Matrix_Type> & aMatrix1,
//               Exp_Temp const & aMatrix2)
//->decltype(aMatrix1.matrix_data()+aMatrix2)
//{
//    return aMatrix1.matrix_data()+aMatrix2;
//}
//
//template<typename Type, typename Matrix_Type >
//auto transpose(Mat<Type,Matrix_Type>  & aMatrix)
//->decltype(aMatrix.matrix_data().transpose())
//{
//    return aMatrix.matrix_data().transpose();
//}
//
//template<typename Type, typename Matrix_Type >
//auto determinant(Mat<Type,Matrix_Type> & aMatrix)
//->decltype(aMatrix.matrix_data().determinant())
//{
//    return aMatrix.matrix_data().determinant();
//}
//
//}
#endif /* SRC_LINALG_CL_XTK_MATRIX_HPP_ */
