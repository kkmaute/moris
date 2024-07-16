/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Bar_3.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_3_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_3_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"    //MTK/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_3 >::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_3 >::get_number_of_points()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_3 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            aIntegrationPoints = {
                { -0.7745966692414834,
                        0.000000000000000,
                        0.7745966692414834 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_3 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 5.555555555555556e-01,
                        8.888888888888888e-01,
                        5.555555555555556e-01 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_3_HPP_ */
