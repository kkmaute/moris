/*
 * cl_MTK_Integration_Coeffs_Quad_4x4.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HEX_2x2x2_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HEX_2x2x2_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp" //LNA/src
#include "cl_MTK_Enums.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::HEX_2x2x2>::get_number_of_dimensions()
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::HEX_2x2x2 >::get_number_of_points()
            {
                return 8;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::HEX_2x2x2 >::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                {
                    -0.577350269189626,
                     0.577350269189626,
                     0.577350269189626,
                    -0.577350269189626,
                    -0.577350269189626,
                     0.577350269189626,
                     0.577350269189626,
                    -0.577350269189626
                 },
                 {
                    -0.577350269189626,
                    -0.577350269189626,
                     0.577350269189626,
                     0.577350269189626,
                    -0.577350269189626,
                    -0.577350269189626,
                     0.577350269189626,
                     0.577350269189626
                 },
                 {
                    -0.577350269189626,
                    -0.577350269189626,
                    -0.577350269189626,
                    -0.577350269189626,
                     0.577350269189626,
                     0.577350269189626,
                     0.577350269189626,
                     0.577350269189626
                 }
            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::HEX_2x2x2 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                    {
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000,
                        1.0000000000000000
                    }
                };
            }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_HEX_2X2X2_HPP_ */