/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Matrix_Base.hpp
 *
 */

#ifndef INCLUDE_CL_LINALG_MATRIX_HPP_
#define INCLUDE_CL_LINALG_MATRIX_HPP_

//#include <memory> // for unique_ptr
//#include <initializer_list>
//
//namespace xtk
//{
//template<typename Matrix_Type>
//class Matrix_Base
//{
//public:
//
//
//    virtual ~Matrix_Base()
//    {
//    }
//
//    //Required Container like properties
//    virtual void resize(size_t aNumRows, size_t aNumColumns) = 0;
//
//    virtual void set_size(size_t aNumRows, size_t aNumColumns) = 0;
//
//    virtual void fill(Type aFillValue) = 0;
//
//    virtual void set_row(size_t aRowIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aRow) = 0;
//
//    virtual void set_column(size_t aColumnIndex, const xtk::Matrix_Base<Type, Matrix_Type> & aColumn) = 0;
//
//    virtual void get_row(size_t aRowIndex, xtk::Matrix_Base<Type, Matrix_Type> & aRow) const = 0;
//
//    virtual void get_column(size_t aColumnIndex, xtk::Matrix_Base<Type, Matrix_Type> & aColumn) const = 0;
//
//    virtual size_t n_rows() const = 0;
//
//    virtual size_t n_cols() const = 0;
//
//    virtual size_t get_num_matrix_elements() const = 0;
//
//    virtual const Type*  data() const = 0;
//
//    virtual Matrix_Type & matrix_data() = 0;
//
//    virtual std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(size_t aNumRows = 0, size_t aNumColumns = 0, Type aFillValue = 0) const = 0;
//
//    virtual std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(std::initializer_list< std::initializer_list< Type > > const & list) const = 0;
//
//    virtual std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> create(Matrix_Type aBaseMatrix) const = 0;
//
//    virtual std::shared_ptr<xtk::Matrix_Base<Type, Matrix_Type>> copy() const  = 0;
//
//    virtual Type get_max_value() const = 0;
//
//    virtual Type get_min_value() const = 0;
//
//    virtual Type & operator()(size_t aRowIndex, size_t aColumnIndex) = 0;
//
//    virtual const Type & operator()(size_t aRowIndex, size_t aColumnIndex) const = 0;
//
//};
//}

#endif /* INCLUDE_CL_LINALG_MATRIX_HPP_ */

