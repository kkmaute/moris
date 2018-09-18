/*
 * cl_XTK_Quad_4_Basis_Function.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_
#define SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_


#include <memory>

#include "../../include/linalg/cl_XTK_Matrix_Base.hpp"
#include "topology/cl_XTK_Basis_Function.hpp"



namespace xtk
{
template<typename Real, typename Real_Matrix>
class Quad_4_Basis_Function: public Basis_Function<Real,Real_Matrix>
{
public:
    Quad_4_Basis_Function()
    {

    }

    void evaluate_basis_function(moris::Matrix< Real_Matrix > const & aLocalCoordinate,
                                 moris::Matrix< Real_Matrix > & aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,4);

        Real tXi = aLocalCoordinate(0,0);
        Real tEta = aLocalCoordinate(0,1);

        aBasisFunctionValues(0,0) = 0.25*(1.0-tXi)*(1.0-tEta);
        aBasisFunctionValues(0,1) = 0.25*(1.0+tXi)*(1.0-tEta);
        aBasisFunctionValues(0,2) = 0.25*(1.0+tXi)*(1.0+tEta);
        aBasisFunctionValues(0,3) = 0.25*(1.0-tXi)*(1.0+tEta);
    }


};
}


#endif /* SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_ */
