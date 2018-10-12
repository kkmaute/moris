/*
 * cl_Matrix_Eigen_Dynamic.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X3_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X3_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"


namespace moris
{

template<typename Type>
class Matrix<Eigen::Matrix<Type, 3, 3>>
{
private:
    Eigen::Matrix<Type,3,3> mMatrix;

public:
    typedef Type Data_Type;

    Matrix(){};

    Matrix(size_t const & aNumRows,
           size_t const & aNumCols)
    {
        MORIS_ASSERT( aNumRows == 3, "Number of rows has to be 3 for 3x3 vector");
        MORIS_ASSERT( aNumCols == 3, "Number of cols has to be 3 for 3x3 vector");
    }

    // template constructor
    template< typename A >
    Matrix(A const & X ):
                mMatrix(X)
     {

     }

    Eigen::Matrix<Type, 3, 3> &
    matrix_data()
    {
        return mMatrix;
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

#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_3X3_HPP_ */
