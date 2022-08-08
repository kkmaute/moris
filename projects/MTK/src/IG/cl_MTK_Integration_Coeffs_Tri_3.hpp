/*
 * cl_MTK_Integration_Coeffs_Tri_3.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_3_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_3_HPP_

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
                Integration_Order::TRI_3 >::get_number_of_dimensions()
        {
            return 2;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3 >::get_number_of_points()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { 0.666666666666667, 0.166666666666667, 0.166666666666667 },
                { 0.166666666666667, 0.666666666666667, 0.166666666666667 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_3 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 0.333333333333333, 0.333333333333333, 0.333333333333333 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_3_HPP_ */
