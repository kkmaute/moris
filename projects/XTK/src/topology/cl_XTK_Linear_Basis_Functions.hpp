/*
 * cl_XTK_Linear_Basis_Functions.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_
#define SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_

namespace xtk
{
template<typename Real, typename Real_Matrix>
class Linear_Basis_Function: public Basis_Function<Real,Real_Matrix>
{
public:
    Linear_Basis_Function()
    {

    }

    void evaluate_basis_function(moris::Matrix< Real_Matrix > const & aLocalCoordinate,
                                 moris::Matrix< Real_Matrix > & aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,2);

        Real tXi = aLocalCoordinate(0,0);

        aBasisFunctionValues(0,0) = 0.5*(1-tXi);
        aBasisFunctionValues(0,1) = 0.5*(1+tXi);
    }


};
}



#endif /* SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_ */
