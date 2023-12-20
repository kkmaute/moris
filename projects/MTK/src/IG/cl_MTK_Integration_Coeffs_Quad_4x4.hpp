/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Quad_4x4.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_4X4_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_4X4_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "moris_typedefs.hpp" //MRS/COR/src
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
            Integration_Order::QUAD_4x4>::get_number_of_dimensions()
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_4x4>::get_number_of_points()
            {
                return 16;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_4x4>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                {  -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053
                 },
                 { -0.861136311594053,
                   -0.861136311594053,
                   -0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                   -0.339981043584856,
                   -0.339981043584856,
                   -0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                    0.861136311594053,
                    0.861136311594053,
                    0.861136311594053
                   }
            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_4x4 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                   { 0.121002993285602,
                     0.226851851851852,
                     0.226851851851852,
                     0.121002993285602,
                     0.226851851851852,
                     0.425293303010694,
                     0.425293303010694,
                     0.226851851851852,
                     0.226851851851852,
                     0.425293303010694,
                     0.425293303010694,
                     0.226851851851852,
                     0.121002993285602,
                     0.226851851851852,
                     0.226851851851852,
                     0.121002993285602
                   }
                };
            }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_4X4_HPP_ */
