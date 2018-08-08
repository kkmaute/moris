/*
 * cl_FEM_Integration_Coeffs_Bar_1.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_5_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_5_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src
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
            Integration_Order::BAR_5>::get_number_of_dimensions()
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_5>::get_number_of_points()
        {
                return 5;
        }

//------------------------------------------------------------------------------

        template<>
        Mat< real >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_5>::get_points()
        {
            Mat< real > aIntegrationPoints =
            {
                {
                    -9.061798459386640e-01,
                    -5.384693101056831e-01,
                     0.000000000000000e+00,
                     5.384693101056831e-01,
                     9.061798459386640e-01
                }
            };

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Mat< real >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_5 >::get_weights()
            {
                Mat< real > aWeights =
                {
                    {
                      2.369268850561891e-01,
                      4.786286704993665e-01,
                      5.688888888888890e-01,
                      4.786286704993665e-01,
                      2.369268850561891e-01
                    }
                };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_5_HPP_ */
