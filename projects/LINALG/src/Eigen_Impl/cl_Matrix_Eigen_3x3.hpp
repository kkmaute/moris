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
class Matrix<Type, Eigen::Matrix<Type, 3, 3>>
{
private:
    Eigen::Matrix<Type,3,3> mMatrix;

public:


    Matrix(size_t const & aNumRows,
            size_t const & aNumCols)
    {

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
