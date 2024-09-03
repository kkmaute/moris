/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Quad_2x2.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_2X2_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_2X2_HPP_

namespace moris::mtk
{
    //------------------------------------------------------------------------------

    template<>
    uint
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_2x2 >::get_number_of_dimensions()
    {
        return 2;
    }

    //------------------------------------------------------------------------------

    template<>
    uint
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_2x2 >::get_number_of_points()
    {
        return 4;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_2x2 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
    {
        aIntegrationPoints = {
            { -0.577350269189626,
                    0.577350269189626,
                    0.577350269189626,
                    -0.577350269189626 },
            { -0.577350269189626,
                    -0.577350269189626,
                    0.577350269189626,
                    0.577350269189626 }
        };
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_2x2 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
    {
        aIntegrationWeights = { { 1.000000000000000,
                1.000000000000000,
                1.000000000000000,
                1.000000000000000 } };
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_2X2_HPP_ */
