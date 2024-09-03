/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Hexahedron_8_Basis_Function.hpp
 *
 */

#ifndef SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_BASIS_FUNCTION_HPP_
#define SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_BASIS_FUNCTION_HPP_

#include "cl_Matrix.hpp"
#include "cl_XTK_Basis_Function.hpp"

namespace moris::xtk
{
    class Hexahedron_8_Basis_Function : public Basis_Function
    {

      public:
        Hexahedron_8_Basis_Function()
        {
        }

        void evaluate_basis_function( Matrix< DDRMat > const &aLocalCoordinate,
                Matrix< DDRMat >                             &aBasisFunctionValues ) const override
        {
            aBasisFunctionValues.resize( 1, 8 );

            moris::real tXi   = aLocalCoordinate( 0 );
            moris::real tEta  = aLocalCoordinate( 1 );
            moris::real tZeta = aLocalCoordinate( 2 );

            aBasisFunctionValues( 0 ) = 0.125 * ( 1.0 - tXi ) * ( 1.0 - tEta ) * ( 1.0 - tZeta );
            aBasisFunctionValues( 1 ) = 0.125 * ( 1.0 + tXi ) * ( 1.0 - tEta ) * ( 1.0 - tZeta );
            aBasisFunctionValues( 2 ) = 0.125 * ( 1.0 + tXi ) * ( 1.0 + tEta ) * ( 1.0 - tZeta );
            aBasisFunctionValues( 3 ) = 0.125 * ( 1.0 - tXi ) * ( 1.0 + tEta ) * ( 1.0 - tZeta );
            aBasisFunctionValues( 4 ) = 0.125 * ( 1.0 - tXi ) * ( 1.0 - tEta ) * ( 1.0 + tZeta );
            aBasisFunctionValues( 5 ) = 0.125 * ( 1.0 + tXi ) * ( 1.0 - tEta ) * ( 1.0 + tZeta );
            aBasisFunctionValues( 6 ) = 0.125 * ( 1.0 + tXi ) * ( 1.0 + tEta ) * ( 1.0 + tZeta );
            aBasisFunctionValues( 7 ) = 0.125 * ( 1.0 - tXi ) * ( 1.0 + tEta ) * ( 1.0 + tZeta );
        }
    };
}    // namespace moris::xtk

#endif /* SRC_TOPOLOGY_CL_XTK_HEXAHEDRON_8_BASIS_FUNCTION_HPP_ */
