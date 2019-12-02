/*
 * cl_XTK_Quad_4_Basis_Function.hpp
 *
 *  Created on: Jul 17, 2017
 *      Author: ktdoble
 */

#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_QUAD_4_BASIS_FUNCTION_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_QUAD_4_BASIS_FUNCTION_HPP_


#include <memory>

#include "../../new/additional/cl_GEN_Basis_Function.hpp"
#include "cl_Matrix.hpp"


namespace moris
{
namespace ge
{
class Quad_4_Basis_Function: public Basis_Function
{
public:
    Quad_4_Basis_Function()
    {

    }

    void evaluate_basis_function(moris::Matrix< moris::DDRMat > const & aLocalCoordinate,
                                 moris::Matrix< moris::DDRMat > & aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,4);

        moris::real tXi = aLocalCoordinate(0,0);
        moris::real tEta = aLocalCoordinate(0,1);

        aBasisFunctionValues(0,0) = 0.25*(1.0-tXi)*(1.0-tEta);
        aBasisFunctionValues(0,1) = 0.25*(1.0+tXi)*(1.0-tEta);
        aBasisFunctionValues(0,2) = 0.25*(1.0+tXi)*(1.0+tEta);
        aBasisFunctionValues(0,3) = 0.25*(1.0-tXi)*(1.0+tEta);
    }


};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_QUAD_4_BASIS_FUNCTION_HPP_ */
