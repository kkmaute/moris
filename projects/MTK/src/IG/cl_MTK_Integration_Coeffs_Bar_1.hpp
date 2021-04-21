/*
 * cl_MTK_Integration_Coeffs_Bar_1.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_1_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_1_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

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
        Integration_Order::BAR_1>::get_number_of_dimensions()
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_1>::get_number_of_points()
        {
                return 1;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs< Integration_Type::GAUSS,
                            Integration_Order::BAR_1 >::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_1 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights.set_size( 1, 1, 2.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_1_HPP_ */