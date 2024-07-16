/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Bar_5.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_5_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_5_HPP_

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
                Integration_Order::BAR_5 >::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_5 >::get_number_of_points()
        {
            return 5;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_5 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            aIntegrationPoints = {
                { -9.061798459386640e-01,
                        -5.384693101056831e-01,
                        0.000000000000000e+00,
                        5.384693101056831e-01,
                        9.061798459386640e-01 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_5 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 2.369268850561891e-01,
                        4.786286704993665e-01,
                        5.688888888888890e-01,
                        4.786286704993665e-01,
                        2.369268850561891e-01 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_5_HPP_ */
