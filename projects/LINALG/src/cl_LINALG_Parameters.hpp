/*
 * cl_LINALG_Parameters.hpp
 *
 *  Created on: Aug 23, 2018
 *      Author: doble
 */

#ifndef PROJECTS_LINALG_SRC_CL_LINALG_PARAMETERS_HPP_
#define PROJECTS_LINALG_SRC_CL_LINALG_PARAMETERS_HPP_


#include "cl_LINALG_Enums.hpp"

namespace moris
{
class Linalg_Parameters
{
public:
    Linalg_Parameters()
    {

    }

    Linalg_Parameters(enum Backend_Dense_Matrix aMatrixBackend):
        mDenseMatrixBackend(aMatrixBackend)
    {

    }

    enum Backend_Dense_Matrix const &
    get_backend_dense_matrix()
    {
        return mDenseMatrixBackend;
    }

private:

    enum Backend_Dense_Matrix mDenseMatrixBackend;
};
}




#endif /* PROJECTS_LINALG_SRC_CL_LINALG_PARAMETERS_HPP_ */
