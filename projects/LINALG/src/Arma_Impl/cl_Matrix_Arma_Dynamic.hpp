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


    Mat_New(size_t const & aNumRows,
            size_t const & aNumCols):
                mMatrix(aNumRows,aNumCols)
    {

    }

    Mat_New(arma::Mat<Type> & aMat):
                mMatrix(aMat)
    {

    }

    // template constructor
    template< typename A >
    Mat_New(
            A const & X ):
            mMatrix(X)
            {

            }


    template< typename ET >
    const Mat_New< Type,  arma::Mat<Type>> &
    operator=(
            ET const & X )
    {
        mMatrix = X;

        return *this;
    }

    arma::Mat<Type> &
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





#endif /* PROJECTS_LINALG_SRC_ARMA_IMPL_CL_MATRIX_ARMA_DYNAMIC_HPP_ */
