/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tri_1.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_1_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_1_HPP_

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
                Integration_Order::TRI_1 >::get_number_of_dimensions()
        {
            return 2;
        }
        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_1 >::get_number_of_points()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_1 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { 0.333333333333333 },
                { 0.333333333333333 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_1 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights.set_size( 1, 1, 1.0 );
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_1_HPP_ */

