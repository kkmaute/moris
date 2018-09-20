/*
 * cl_FEM_Integration_Coeffs_Quad_4x4.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_HEX_2x2x2_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_HEX_2x2x2_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
#include "linalg_typedefs.hpp" //LNA/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
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
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::HEX_2x2x2 >::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints =
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

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::HEX_2x2x2 >::get_weights()
            {
                Matrix< DDRMat > aWeights =
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

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_HEX_2X2X2_HPP_ */
