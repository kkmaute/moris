/*
 * cl_XTK_Parameters.hpp
 *
 *  Created on: Jun 15, 2018
 *      Author: ktdoble
 */

#ifndef SRC_CORE_CL_XTK_PARAMETERS_HPP_
#define SRC_CORE_CL_XTK_PARAMETERS_HPP_


#include "xtk/cl_XTK_Enums.hpp"
namespace xtk
{
class Parameters
{
public:
    Parameters()
    {

    }

    Parameters(enum Matrix_Backend aMatrixBackend):
        mMatrixBackend(aMatrixBackend)
    {

    }


    enum Matrix_Backend const &
    get_matrix_backend()
    {
        return mMatrixBackend;
    }

private:

    enum Matrix_Backend mMatrixBackend;
};
}



#endif /* SRC_CORE_CL_XTK_PARAMETERS_HPP_ */
