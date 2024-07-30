/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Bar_6.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_16_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_16_HPP_

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
                Integration_Order::BAR_16 >::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_16 >::get_number_of_points()
        {
            return 16;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_16 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            aIntegrationPoints = {
                { -0.9894009349916499,
                        -0.9445750230732326,
                        -0.8656312023878318,
                        -0.7554044083550030,
                        -0.6178762444026438,
                        -0.4580167776572274,
                        -0.2816035507792589,
                        -0.0950125098376374,
                        0.0950125098376374,
                        0.2816035507792589,
                        0.4580167776572274,
                        0.6178762444026438,
                        0.7554044083550030,
                        0.8656312023878318,
                        0.9445750230732326,
                        0.9894009349916499 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_16 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 0.0271524594117541,
                        0.0622535239386479,
                        0.0951585116824928,
                        0.1246289712555339,
                        0.1495959888165767,
                        0.1691565193950025,
                        0.1826034150449236,
                        0.1894506104550685,
                        0.1894506104550685,
                        0.1826034150449236,
                        0.1691565193950025,
                        0.1495959888165767,
                        0.1246289712555339,
                        0.0951585116824928,
                        0.0622535239386479,
                        0.0271524594117541 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_16_HPP_ */
