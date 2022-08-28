/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Quad_4_Basis_Function.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_
#define SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_

#include <memory>

#include "cl_Matrix.hpp"
#include "cl_XTK_Basis_Function.hpp"

namespace xtk
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

            moris::real tXi = aLocalCoordinate(0);
            moris::real tEta = aLocalCoordinate(1);

            aBasisFunctionValues(0) = 0.25*(1.0-tXi)*(1.0-tEta);
            aBasisFunctionValues(1) = 0.25*(1.0+tXi)*(1.0-tEta);
            aBasisFunctionValues(2) = 0.25*(1.0+tXi)*(1.0+tEta);
            aBasisFunctionValues(3) = 0.25*(1.0-tXi)*(1.0+tEta);
        }
    };
}

#endif /* SRC_TOPOLOGY_CL_XTK_QUAD_4_BASIS_FUNCTION_HPP_ */

