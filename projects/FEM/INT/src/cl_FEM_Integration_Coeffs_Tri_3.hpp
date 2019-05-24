/*
 * cl_FEM_Integration_Coeffs_Tri_3.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_

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
            Integration_Order::TRI_3>::get_number_of_dimensions()
        {
            return 3;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_3>::get_number_of_points()
            {
                return 3;
            }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3>::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints =
            {
                { 0.666666666666667, 0.166666666666667, 0.166666666666667 },
                { 0.166666666666667, 0.666666666666667, 0.166666666666667 },
                { 0.166666666666667, 0.166666666666667, 0.666666666666667 }
            };

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3 >::get_weights()
            {
                Matrix< DDRMat > aWeights =
                {
                    { 0.333333333333333, 0.333333333333333,0.333333333333333 }
                };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_ */
