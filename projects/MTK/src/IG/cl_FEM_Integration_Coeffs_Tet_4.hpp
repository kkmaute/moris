/*
 * cl_FEM_Integration_Coeffs_Tet_4.hpp
 *
 *  Created on: Apr 05, 2019
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_4_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_4_HPP_

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
            Integration_Order::TET_4>::get_number_of_dimensions()
        {
            return 4;
        }
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TET_4>::get_number_of_points()
            {
                return 4;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_4>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                {
                    0.585410196624968,
                    0.138196601125010,
                    0.138196601125010,
                    0.138196601125010
                },
                {
                    0.138196601125010,
                    0.585410196624968,
                    0.138196601125010,
                    0.138196601125010
                },
                {
                    0.138196601125010,
                    0.138196601125010,
                    0.585410196624968,
                    0.138196601125010
                },
                {
                    0.138196601125010,
                    0.138196601125010,
                    0.138196601125010,
                    0.585410196624968
                }
            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_4 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                    {
                        0.250000000000000,
                        0.250000000000000,
                        0.250000000000000,
                        0.250000000000000
                    }
                };
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TET_4_HPP_ */
