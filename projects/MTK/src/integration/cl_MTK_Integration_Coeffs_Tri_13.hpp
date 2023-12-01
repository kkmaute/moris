/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tri_13.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_13_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_13_HPP_

// MRS/COR/src
#include "typedefs.hpp"
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
                Integration_Order::TRI_13 >::get_number_of_dimensions()
        {
            return 2;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_13 >::get_number_of_points()
        {
            return 13;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_13 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { 0.333333333333333,

                        0.479308067841920,
                        0.260345966079040,
                        0.260345966079040,

                        0.869739794195568,
                        0.065130102902216,
                        0.065130102902216,

                        0.048690315425316,
                        0.048690315425316,
                        0.312865496004874,
                        0.312865496004874,
                        0.638444188569810,
                        0.638444188569810 },
                { 0.333333333333333,

                        0.260345966079040,
                        0.479308067841920,
                        0.260345966079040,

                        0.065130102902216,
                        0.869739794195568,
                        0.065130102902216,

                        0.312865496004874,
                        0.638444188569810,
                        0.048690315425316,
                        0.638444188569810,
                        0.048690315425316,
                        0.312865496004874 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_13 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { -0.149570044467682,

                        0.175615257433208,
                        0.175615257433208,
                        0.175615257433208,

                        0.053347235608838,
                        0.053347235608838,
                        0.053347235608838,

                        0.077113760890257,
                        0.077113760890257,
                        0.077113760890257,
                        0.077113760890257,
                        0.077113760890257,
                        0.077113760890257 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_13_HPP_ */

