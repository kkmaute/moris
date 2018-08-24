/*
 * cl_LINALG_Backend_Factory.hpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_CL_DENSE_BACKEND_MATRIX_FACTORY_HPP_
#define PROJECTS_LINALG_SRC_CL_DENSE_BACKEND_MATRIX_FACTORY_HPP_

#include <memory>
#include "cl_Matrix_Base.hpp" //LINALG/src
#include "cl_LINALG_Parameters.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/cl_Matrix_Base_Eigen.hpp" //LINALG/src
#endif
#ifdef MORIS_USE_ARMA
#include "Arma_Impl/cl_Matrix_Base_Arma.hpp" //LINALG/src

namespace moris
{

template<typename T1, typename T2>
class A {};            // primary template

template<typename Type>
class A<Type, arma::Mat<Type>>
{
public:
        std::shared_ptr<Matrix_Base<Type,arma::Mat<Type>>>
        create_matrix_base(size_t const & aNumRows,
                           size_t const & aNumCols)
        {
            return std::make_shared<Matrix_Base_Arma_Dynamic<Type>>(aNumRows,aNumCols);
        }
};

#endif

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/cl_Matrix_Base_Eigen.hpp" //LINALG/src

template<typename Type>
class A<Type, Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic>>
{
public:
    std::shared_ptr<Matrix_Base<Type,Eigen::Matrix<Type,Eigen::Dynamic,Eigen::Dynamic>>>
    create_matrix_base(size_t const & aNumRows,
                       size_t const & aNumCols)
    {
        return std::make_shared<Matrix_Base_Eigen_Dynamic<Type>>(aNumRows,aNumCols);
    }
};
#endif
}

#endif /* PROJECTS_LINALG_SRC_CL_DENSE_BACKEND_MATRIX_FACTORY_HPP_ */
