/*
 * cl_Matrix_Eigen_Dynamic.hpp
 *
 *  Created on: Aug 24, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_
#define PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_

#include "cl_Matrix.hpp"
#include "Eigen/Dense"

namespace moris
{
template<typename Type>
class Mat_New<Type, Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>>
{
private:
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> mMatrix;

public:
    Mat_New()
    {

    };

    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols):
                mMatrix(aNumRows,aNumCols)
    {

    }
    // template constructor
    template< typename A >
    Mat_New(A const & X ):
                mMatrix(X)
     {

     }

    inline
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> &
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


#endif /* PROJECTS_LINALG_SRC_EIGEN_IMPL_CL_MATRIX_EIGEN_DYNAMIC_HPP_ */
