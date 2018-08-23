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


namespace moris
{

/*
 * This class wraps the matrix_base interface and hides the shared pointer.
 * All functions added to linalg/cl_XTK_Matrix_Base.hpp must be wrapped here
 */

// Global parameters
extern Parameters gParameters;

template < typename Type,
           typename Matrix_Type >
class Dense_Backend_Matrix_Factory
{

public:
    Dense_Backend_Matrix_Factory(){};

    std::shared_ptr<Matrix_Base<Type,Matrix_Type>>
    create_matrix_base(size_t const & aNumRows,
                       size_t const & aNumCols )
    {
        enum Matrix_Backend tMatrixBackend = gParameters.get_matrix_backend();

        std::shared_ptr<Matrix_Base<Type,Matrix_Type>> tMatrixBase;
        switch(tMatrixBackend)
        {
#ifdef MORIS_USE_EIGEN
        case(Backend_Dense_Matrix::EIGEN):
        {
            tMatrixBase = std::make_shared<Matrix_Base_Eigen<Type,Matrix_Type>>(aNumRows,aNumCols);
            break;
        }
#endif

#ifdef MORIS_USE_ARMA
#endif
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

#endif /* PROJECTS_LINALG_SRC_CL_DENSE_BACKEND_MATRIX_FACTORY_HPP_ */
