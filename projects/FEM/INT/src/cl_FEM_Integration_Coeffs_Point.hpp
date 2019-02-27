/*
 * cl_FEM_Integration_Coeffs_Point.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
        Integration_Order::POINT>::get_number_of_dimensions()
        {
            return 0;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::POINT>::get_number_of_points()
        {
                return 1;
        }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Integration_Coeffs< Integration_Type::GAUSS,
                            Integration_Order::POINT >::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints =
            {
                { 0.0 }
            };

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::POINT >::get_weights()
            {
                Matrix< DDRMat > aWeights =
                {
                    { 1.0 }
                };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_ */
