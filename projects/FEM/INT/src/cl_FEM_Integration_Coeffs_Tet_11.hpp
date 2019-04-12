/*
 * cl_FEM_Integration_Coeffs_Tet_11.hpp
 *
 *  Created on: Apr 05, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_11_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_11_HPP_

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
            Integration_Order::TET_11>::get_number_of_dimensions()
        {
            return 4;
        }
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TET_11>::get_number_of_points()
            {
                return 11;
            }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11>::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints =
            {
                { 0.7857142857142860, 0.0714285714285714, 0.0714285714285714, 0.0714285714285714,
                  0.3994035761667990, 0.3994035761667990, 0.3994035761667990, 0.1005964238332010,
                  0.1005964238332010, 0.1005964238332010, 0.250000000000000 },
                { 0.0714285714285714, 0.7857142857142860, 0.0714285714285714, 0.0714285714285714,
                  0.3994035761667990, 0.1005964238332010, 0.1005964238332010, 0.3994035761667990,
                  0.3994035761667990, 0.1005964238332010, 0.250000000000000 },
                { 0.0714285714285714, 0.0714285714285714, 0.7857142857142860, 0.0714285714285714,
                  0.1005964238332010, 0.3994035761667990, 0.1005964238332010, 0.3994035761667990,
                  0.1005964238332010, 0.3994035761667990, 0.250000000000000 },
                { 0.0714285714285714, 0.0714285714285714, 0.0714285714285714, 0.7857142857142860,
                  0.1005964238332010, 0.1005964238332010, 0.3994035761667990, 0.1005964238332010,
                  0.3994035761667990, 0.3994035761667990, 0.250000000000000 }
            };
            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11>::get_weights()
            {
                Matrix< DDRMat > aWeights =
                {
                    { 0.0076222222222222, 0.0076222222222222, 0.0076222222222222, 0.0076222222222222,
                      0.0248888888888889, 0.0248888888888889, 0.0248888888888889, 0.0248888888888889,
                      0.0248888888888889, 0.0248888888888889,-0.013155555555555 }
                };
                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_11_HPP_ */
