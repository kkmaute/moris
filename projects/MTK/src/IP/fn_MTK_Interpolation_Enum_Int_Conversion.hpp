/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_MTK_Interpolation_Enum_Int_Conversion.hpp
 *
 */

#ifndef SRC_fn_MTK_Interpolation_Enum_Int_Conversion
#define SRC_fn_MTK_Interpolation_Enum_Int_Conversion

#include "assert.hpp"
#include "typedefs.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris
{
    namespace mtk
    {
        //------------------------------------------------------------------------------

        inline uint
        ip_order_enum_to_uint( enum Interpolation_Order aIpOrder )
        {
            switch ( aIpOrder )
            {
                case Interpolation_Order::CONSTANT:
                {
                    return 0;
                    break;
                }

                case Interpolation_Order::LINEAR:
                {
                    return 1;
                    break;
                }

                case Interpolation_Order::QUADRATIC:
                {
                    return 2;
                    break;
                }

                case Interpolation_Order::CUBIC:
                {
                    return 3;
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "ip_order_enum_to_uint() - aIpOrder doesn't define polynomial order / unknown." );
                    return 0;
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

        inline enum Interpolation_Order
        ip_order_uint_to_enum( uint aIpOrder )
        {
            switch ( aIpOrder )
            {
                case 0:
                {
                    return Interpolation_Order::CONSTANT;
                    break;
                }

                case 1:
                {
                    return Interpolation_Order::LINEAR;
                    break;
                }

                case 2:
                {
                    return Interpolation_Order::QUADRATIC;
                    break;
                }

                case 3:
                {
                    return Interpolation_Order::CUBIC;
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "ip_order_uint_to_enum() - Only interpolation orders 0-3 supported." );
                    return Interpolation_Order::UNDEFINED;
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace mtk
}    // namespace moris

#endif /* fn_MTK_Interpolation_Enum_Int_Conversion.hpp */
