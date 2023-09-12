/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Hexahedron_8_Basis_Function.hpp
 *
 */

#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_HEXAHEDRON_8_BASIS_FUNCTION_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_HEXAHEDRON_8_BASIS_FUNCTION_HPP_

#include "../../new/additional/cl_GEN_Basis_Function.hpp"
#include "cl_Matrix.hpp"

namespace moris
{
namespace ge
{
class Hexahedron_8_Basis_Function: public Basis_Function
{
public:
    Hexahedron_8_Basis_Function()
    {

    }

    void evaluate_basis_function(moris::Matrix< moris::DDRMat > const & aLocalCoordinate,
                                 moris::Matrix< moris::DDRMat > & aBasisFunctionValues) const
    {
        aBasisFunctionValues.resize(1,8);

        moris::real tXi = aLocalCoordinate(0,0);
        moris::real tEta = aLocalCoordinate(0,1);
        moris::real tZeta = aLocalCoordinate(0,2);

        aBasisFunctionValues(0,0) = 0.125*(1.0-tXi)*(1.0-tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,1) = 0.125*(1.0+tXi)*(1.0-tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,2) = 0.125*(1.0+tXi)*(1.0+tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,3) = 0.125*(1.0-tXi)*(1.0+tEta)*(1.0-tZeta);
        aBasisFunctionValues(0,4) = 0.125*(1.0-tXi)*(1.0-tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,5) = 0.125*(1.0+tXi)*(1.0-tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,6) = 0.125*(1.0+tXi)*(1.0+tEta)*(1.0+tZeta);
        aBasisFunctionValues(0,7) = 0.125*(1.0-tXi)*(1.0+tEta)*(1.0+tZeta);
    }

};
}
}

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_HEXAHEDRON_8_BASIS_FUNCTION_HPP_ */

