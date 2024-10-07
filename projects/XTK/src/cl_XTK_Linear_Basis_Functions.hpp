/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Linear_Basis_Functions.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_
#define SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_

namespace moris::xtk
{
    class Linear_Basis_Function : public Basis_Function
    {

      public:
        Linear_Basis_Function()
        {
        }

        void evaluate_basis_function( Matrix< DDRMat > const &aLocalCoordinate,
                Matrix< DDRMat >                             &aBasisFunctionValues ) const override
        {
            aBasisFunctionValues.resize( 1, 2 );

            moris::real tXi = aLocalCoordinate( 0 );

            aBasisFunctionValues( 0 ) = 0.5 * ( 1 - tXi );
            aBasisFunctionValues( 1 ) = 0.5 * ( 1 + tXi );
        }
    };
}    // namespace moris::xtk

#endif /* SRC_TOPOLOGY_CL_XTK_LINEAR_BASIS_FUNCTIONS_HPP_ */
