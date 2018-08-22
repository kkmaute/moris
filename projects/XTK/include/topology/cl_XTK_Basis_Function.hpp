/*
 * cl_XTK_Basis_Function.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef INCLUDE_TOPOLOGY_CL_XTK_BASIS_FUNCTION_HPP_
#define INCLUDE_TOPOLOGY_CL_XTK_BASIS_FUNCTION_HPP_

#include "linalg/cl_XTK_Matrix.hpp"

namespace xtk
{
template<typename Real, typename Real_Matrix>
class Basis_Function
{
public:
    ~Basis_Function()
    {

    }

    /**
     * Evaluates the values of the basis functions at a give local or parametric coordinate [-1,1]
     */
    virtual void evaluate_basis_function(Mat<Real,Real_Matrix> const & aLocalCoordinate,
                                         Mat<Real,Real_Matrix> & aBasisFunctionValues) const = 0;
};

}


#endif /* INCLUDE_TOPOLOGY_CL_XTK_BASIS_FUNCTION_HPP_ */
