/*
 * cl_FEM_Integration_Coeffs_Tet_5.hpp
 *
 *  Created on: Apr 05, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_5_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_5_HPP_

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
            Integration_Order::TET_5>::get_number_of_dimensions()
        {
            return 4;
        }
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TET_5>::get_number_of_points()
            {
                return 5;
            }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5>::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints =
            {
                {0.500000000000000, 0.166666666666667, 0.166666666666667, 0.166666666666667, 0.250000000000000},
                {0.166666666666667, 0.500000000000000, 0.166666666666667, 0.166666666666667, 0.250000000000000},
                {0.166666666666667, 0.166666666666667, 0.500000000000000, 0.166666666666667, 0.250000000000000},
                {0.166666666666667, 0.166666666666667, 0.166666666666667, 0.500000000000000, 0.250000000000000}

            };
            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5 >::get_weights()
            {
                Matrix< DDRMat > aWeights =
                {
                    {0.450000000000000, 0.450000000000000, 0.450000000000000, 0.450000000000000,-0.800000000000000}
                };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_5_HPP_ */
