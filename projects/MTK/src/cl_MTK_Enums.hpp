/*
 * cl_MTK_Enums.hpp
 *
 *  Created on: Jul 25, 2018
 *      Author: messe
 */

#ifndef SRC_MESH_CL_MTK_ENUMS_HPP_
#define SRC_MESH_CL_MTK_ENUMS_HPP_

#include "assert.hpp"
#include "typedefs.hpp"

//------------------------------------------------------------------------------
namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        enum class Geometry_Type
        {
            LINE,  // 1D line or curve
            QUAD,  // rectangle
            TRI,   // triangle
            HEX,   // quadrangle
            TET,   // tetrahedron
            PENTA,  // pentahedron
            UNDEFINED
        };

//------------------------------------------------------------------------------

        enum class Interpolation_Order
        {
            CONSTANT,
            LINEAR,
            QUADRATIC,
            SERENDIPITY,
            CUBIC,
            UNDEFINED
        };

//------------------------------------------------------------------------------

        /**
         * converts an interpolation order to a numeric value
         */
        uint
        interpolation_order_to_uint( const Interpolation_Order& aOrder )
        {
            switch( aOrder )
            {
                case( Interpolation_Order::CONSTANT ) :
                {
                    return 0;
                    break;
                }
                case( Interpolation_Order::LINEAR ) :
                {
                    return 1;
                    break;
                }
                case( Interpolation_Order::SERENDIPITY ) :
                {
                    MORIS_ERROR( false, "Interpolation_Order::SERENDIPITY cannot be converted to uint" );
                    return 0;
                    break;
                }
                case( Interpolation_Order::QUADRATIC ) :
                {
                    return 2;
                    break;

                }
                case( Interpolation_Order::CUBIC ) :
                {
                    return 3;
                    break;
                }
                default :
                {
                    MORIS_ERROR( false, "unknown interpolation order" );
                    return 0;
                    break;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


#endif /* SRC_MESH_CL_MTK_ENUMS_HPP_ */
