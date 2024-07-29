/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Bar_6.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_6_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_6_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Coeffs.hpp"

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_6 >::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_6 >::get_number_of_points()
        {
            return 6;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_6 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            aIntegrationPoints = {
                { -0.932469514203152,
                        -0.661209386466265,
                        -0.238619186083197,
                        0.238619186083197,
                        0.661209386466265,
                        0.932469514203152 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_6 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 0.171324492379170,
                        0.360761573048139,
                        0.467913934572691,
                        0.467913934572691,
                        0.360761573048139,
                        0.171324492379170 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_6_HPP_ */
