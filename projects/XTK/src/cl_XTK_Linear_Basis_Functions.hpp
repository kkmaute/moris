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
class Linear_Basis_Function: public Basis_Function
{
public:
    Linear_Basis_Function()
    {

    }

    void evaluate_basis_function(moris::Matrix< moris::DDRMat > const & aLocalCoordinate,
                                 moris::Matrix< moris::DDRMat > & aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,2);

        moris::real tXi = aLocalCoordinate(0,0);

        aBasisFunctionValues(0,0) = 0.5*(1-tXi);
        aBasisFunctionValues(0,1) = 0.5*(1+tXi);
    }


};
}



#endif /* SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_ */
