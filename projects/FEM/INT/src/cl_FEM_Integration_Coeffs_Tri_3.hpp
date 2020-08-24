/*
 * cl_FEM_Integration_Coeffs_Tri_3.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Integration_Coeffs.hpp"

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
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                { 0.666666666666667, 0.166666666666667, 0.166666666666667 },
                { 0.166666666666667, 0.666666666666667, 0.166666666666667 },
                { 0.166666666666667, 0.166666666666667, 0.666666666666667 }
            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                    { 0.333333333333333, 0.333333333333333,0.333333333333333 }
                };
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_3_HPP_ */
