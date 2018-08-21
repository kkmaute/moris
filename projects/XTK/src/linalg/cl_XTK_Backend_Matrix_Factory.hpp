/*
 * cl_XTK_Backend_Matrix_Factory.hpp
 *
 *  Created on: Aug 2, 2018
 *      Author: ktdoble
 */

#ifndef SRC_LINALG_CL_XTK_BACKEND_MATRIX_FACTORY_HPP_
#define SRC_LINALG_CL_XTK_BACKEND_MATRIX_FACTORY_HPP_

// Default type defs
#include<iostream>

#include "linalg/cl_XTK_Matrix_Base.hpp"

#ifdef XTK_USE_EIGEN
#include "linalg/cl_XTK_Matrix_Eigen.hpp"
#else
#endif


namespace xtk
{

/*
 * This class wraps the matrix_base interface and hides the shared pointer.
 * All functions added to linalg/cl_XTK_Matrix_Base.hpp must be wrapped here
 */

// Global parameters
extern Parameters gParameters;

template < typename Type,
           typename Matrix_Type >
class Backend_Matrix_Factory
{

public:
    Backend_Matrix_Factory(){};

    std::shared_ptr<Matrix_Base<Type,Matrix_Type>>
    create_matrix_base(size_t const & aNumRows,
                       size_t const & aNumCols )
    {
        enum Matrix_Backend tMatrixBackend = gParameters.get_matrix_backend();

        std::shared_ptr<Matrix_Base<Type,Matrix_Type>> tMatrixBase;
        switch(tMatrixBackend)
        {

        case(Matrix_Backend::EIGEN):
        {
#ifdef XTK_USE_EIGEN
            tMatrixBase = std::make_shared<Matrix_Base_Eigen<Type,Matrix_Type>>(aNumRows,aNumCols);
#endif
            break;
        }

        case(Matrix_Backend::TEUCHOS):
        {
#ifndef XTK_USE_EIGEN
            std::cout<<"Teuchos not implemented"<<std::endl;
            break;
#endif
        }

        default:
        {
            std::cout<<"Matrix backend not implemented"<<std::endl;
            break;
        }
        }
        return tMatrixBase;
    }

    std::shared_ptr<Matrix_Base<Type,Matrix_Type>>
    create_matrix_base(std::initializer_list<std::initializer_list<Type> > const & aInitList )
    {
        enum Matrix_Backend tMatrixBackend = gParameters.get_matrix_backend();

        std::shared_ptr<Matrix_Base<Type,Matrix_Type>> tMatrixBase;
        switch(tMatrixBackend)
        {

        case(Matrix_Backend::EIGEN):
        {
#ifdef XTK_USE_EIGEN
            tMatrixBase = std::make_shared<Matrix_Base_Eigen<Type,Matrix_Type>>(aInitList);
#endif
            break;

        }

        case(Matrix_Backend::TEUCHOS):
        {
            std::cout<<"Teuchos not implemented"<<std::endl;
            break;
        }

        default:
        {
            std::cout<<"Matrix backend not implemented"<<std::endl;
            break;
        }
        }
        return tMatrixBase;
    }

    std::shared_ptr<Matrix_Base<Type,Matrix_Type>>
    create_matrix_base(Matrix_Type const & aBackendMat )
    {
        enum Matrix_Backend tMatrixBackend = gParameters.get_matrix_backend();

        std::shared_ptr<Matrix_Base<Type,Matrix_Type>> tMatrixBase;
        switch(tMatrixBackend)
        {

        case(Matrix_Backend::EIGEN):
        {
#ifdef XTK_USE_EIGEN
            tMatrixBase = std::make_shared<Matrix_Base_Eigen<Type,Matrix_Type>>(aBackendMat);
#endif

            break;
        }

        case(Matrix_Backend::TEUCHOS):
        {
            std::cout<<"Teuchos not implemented"<<std::endl;
            break;
        }

        default:
        {
            std::cout<<"Matrix backend not implemented"<<std::endl;
            break;
        }
        }
        return tMatrixBase;
    }

};
}


#endif /* SRC_LINALG_CL_XTK_BACKEND_MATRIX_FACTORY_HPP_ */
