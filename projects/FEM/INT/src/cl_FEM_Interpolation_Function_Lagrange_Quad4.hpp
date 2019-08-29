/*
 * cl_FEM_Interpolation_Function_Quad4.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_

#include "assert.h"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Interpolation_Function.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_number_of_param_dimensions() const
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::LINEAR;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat.set_size( 2, 4, 0.0 );
            aXiHat( 0, 0 ) = -1.000000;
            aXiHat( 1, 0 ) = -1.000000;
            aXiHat( 0, 1 ) =  1.000000;
            aXiHat( 1, 1 ) = -1.000000;
            aXiHat( 0, 2 ) =  1.000000;
            aXiHat( 1, 2 ) =  1.000000;
            aXiHat( 0, 3 ) = -1.000000;
            aXiHat( 1, 3 ) =  1.000000;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                              Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_N: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // populate matrix with values
            aNXi.set_size( 1, 4 );
            aNXi( 0 ) = ( ( 1.0 - xi ) * ( 1.0 - eta ) ) * 0.25;
            aNXi( 1 ) = ( ( 1.0 + xi ) * ( 1.0 - eta ) ) * 0.25;
            aNXi( 2 ) = ( ( 1.0 + xi ) * ( 1.0 + eta ) ) * 0.25;
            aNXi( 3 ) = ( ( 1.0 - xi ) * ( 1.0 + eta ) ) * 0.25;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                  Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_dNdXi: aXi not allocated or hat wrong size." );

            // unpack xi and eta from input vector
            real  xi = aXi( 0 );
            real eta = aXi( 1 );

            // populate output matrix
            adNdXi.set_size( 2, 4 );
            adNdXi( 0, 0 ) =  0.25 * ( eta - 1.0 );
            adNdXi( 1, 0 ) =  0.25 * ( xi - 1.0 );

            adNdXi( 0, 1 ) = -0.25 * ( eta - 1.0 );
            adNdXi( 1, 1 ) = -0.25 * ( xi + 1.0 );

            adNdXi( 0, 2 ) =  0.25 * ( eta + 1.0 );
            adNdXi( 1, 2 ) =  0.25 * ( xi + 1.0 );

            adNdXi( 0, 3 ) = -0.25 * ( eta + 1.0 );
            adNdXi( 1, 3 ) = -0.25 * ( xi - 1.0 );
        }


//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad2NdXi2 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad2NdXi2.set_size( 3, 4 );
            ad2NdXi2( 0, 0 ) =  0.0;
            ad2NdXi2( 1, 0 ) =  0.0;
            ad2NdXi2( 2, 0 ) =  0.25;

            ad2NdXi2( 0, 1 ) =  0.0;
            ad2NdXi2( 1, 1 ) =  0.0;
            ad2NdXi2( 2, 1 ) = -0.25;

            ad2NdXi2( 0, 2 ) =  0.0;
            ad2NdXi2( 1, 2 ) =  0.0;
            ad2NdXi2( 2, 2 ) =  0.25;

            ad2NdXi2( 0, 3 ) =  0.0;
            ad2NdXi2( 1, 3 ) =  0.0;
            ad2NdXi2( 2, 3 ) = -0.25;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::QUAD, Interpolation_Type::LAGRANGE, 2, 4 >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                    Matrix< DDRMat > & ad3NdXi3 ) const
        {

            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 2, "QUAD4 - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            // populate output matrix
            ad3NdXi3.set_size( 4, 4, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_QUAD4_HPP_ */
