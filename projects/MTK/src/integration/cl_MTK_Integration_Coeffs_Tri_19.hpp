/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tri_19.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_19_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_19_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
// LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integration_Coeffs.hpp"

namespace moris::mtk
{
    //------------------------------------------------------------------------------
    template<>
    uint
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_19 >::get_number_of_dimensions()
    {
        return 2;
    }

    //------------------------------------------------------------------------------

    template<>
    uint
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_19 >::get_number_of_points()
    {
        return 19;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_19 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
    {
        aIntegrationPoints = {
            { 0.333333333333333,

                    0.020634961602525,
                    0.489682519198738,
                    0.489682519198738,

                    0.125820817014127,
                    0.437089591492937,
                    0.437089591492937,

                    0.623592928761935,
                    0.188203535619033,
                    0.188203535619033,

                    0.910540973211095,
                    0.044729513394453,
                    0.044729513394453,

                    0.036838412054736,
                    0.036838412054736,
                    0.221962989160766,
                    0.221962989160766,
                    0.741198598784498,
                    0.741198598784498 },
            { 0.333333333333333,

                    0.489682519198738,
                    0.020634961602525,
                    0.489682519198738,

                    0.437089591492937,
                    0.125820817014127,
                    0.437089591492937,

                    0.188203535619033,
                    0.623592928761935,
                    0.188203535619033,

                    0.044729513394453,
                    0.910540973211095,
                    0.044729513394453,

                    0.221962989160766,
                    0.741198598784498,
                    0.036838412054736,
                    0.741198598784498,
                    0.221962989160766,
                    0.036838412054736 }
        };
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_19 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
    {
        aIntegrationWeights = {
            { 0.097135796282799,

                    0.031334700227139,
                    0.031334700227139,
                    0.031334700227139,

                    0.077827541004774,
                    0.077827541004774,
                    0.077827541004774,

                    0.079647738927210,
                    0.079647738927210,
                    0.079647738927210,

                    0.025577675658698,
                    0.025577675658698,
                    0.025577675658698,

                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289,
                    0.043283539377289 }
        };
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_19_HPP_ */
