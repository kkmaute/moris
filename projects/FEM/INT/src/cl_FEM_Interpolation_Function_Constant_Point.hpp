/*
 * cl_FEM_Interpolation_Function_Constant_Point.hpp
 *
 *  Created on: May 30, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_CONSTANT_POINT_HPP_
#define SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_CONSTANT_POINT_HPP_

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
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1 >::get_number_of_param_dimensions() const
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        mtk::Interpolation_Order
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1 >::get_interpolation_order() const
        {
            return mtk::Interpolation_Order::UNDEFINED;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1 >::get_param_coords( Matrix< DDRMat > & aXiHat ) const
        {
            aXiHat.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1  >::eval_N( const Matrix< DDRMat > & aXi,
                                                                                                                Matrix< DDRMat > & aNXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "POINT - eval_N: aXi not allocated or hat wrong size." );

            aNXi.set_size( 1, 1, 1.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1 >::eval_dNdXi( const Matrix< DDRMat > & aXi,
                                                                                                                   Matrix< DDRMat > & adNdXi ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "POINT - eval_dNdXi: aXi not allocated or hat wrong size." );

            adNdXi.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1  >::eval_d2NdXi2( const Matrix< DDRMat > & aXi,
                                                                                                                      Matrix< DDRMat > & ad2NdXi2 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "POINT - eval_d2NdXi2: aXi not allocated or hat wrong size." );

            ad2NdXi2.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Interpolation_Function< mtk::Geometry_Type::POINT, Interpolation_Type::CONSTANT, 1, 1  >::eval_d3NdXi3( const Matrix< DDRMat > & aXi,
                                                                                                                      Matrix< DDRMat > & ad3NdXi3 ) const
        {
            // make sure that input is correct
            MORIS_ASSERT( aXi.length() >= 1, "POINT - eval_d3NdXi3: aXi not allocated or hat wrong size." );

            ad3NdXi3.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

//------------------------------------------------------------------------------
#endif /* SRC_FEM_CL_FEM_INTERPOLATION_FUNCTION_LAGRANGE_BAR1_HPP_ */
