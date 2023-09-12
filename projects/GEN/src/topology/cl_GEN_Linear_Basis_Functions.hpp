/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Linear_Basis_Functions.hpp
 *
 */

#ifndef PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_LINEAR_BASIS_FUNCTIONS_HPP_
#define PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_LINEAR_BASIS_FUNCTIONS_HPP_

namespace moris
{
namespace ge
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
}

#endif /* PROJECTS_GEN_SRC_NEW_ADDITIONAL_CL_GEN_LINEAR_BASIS_FUNCTIONS_HPP_ */

