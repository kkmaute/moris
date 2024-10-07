/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tri_12.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_12_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_12_HPP_

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
            Integration_Order::TRI_12 >::get_number_of_dimensions()
    {
        return 2;
    }
    //------------------------------------------------------------------------------

    template<>
    uint
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_12 >::get_number_of_points()
    {
        return 12;
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_12 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
    {
        aIntegrationPoints = {
            { 0.873821971016996,
                    0.063089014491502,
                    0.063089014491502,

                    0.501426509658179,
                    0.249286745170910,
                    0.249286745170910,

                    0.636502499121399,
                    0.636502499121399,
                    0.310352451033785,
                    0.310352451033785,
                    0.053145049844816,
                    0.053145049844816 },
            { 0.063089014491502,
                    0.873821971016996,
                    0.063089014491502,

                    0.249286745170910,
                    0.501426509658179,
                    0.249286745170910,

                    0.310352451033785,
                    0.053145049844816,
                    0.636502499121399,
                    0.053145049844816,
                    0.310352451033785,
                    0.636502499121399 }
        };
    }

    //------------------------------------------------------------------------------

    template<>
    void
    Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_12 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
    {
        aIntegrationWeights = {
            { 0.050844906370207,
                    0.050844906370207,
                    0.050844906370207,

                    0.116786275726379,
                    0.116786275726379,
                    0.116786275726379,

                    0.082851075618374,
                    0.082851075618374,
                    0.082851075618374,
                    0.082851075618374,
                    0.082851075618374,
                    0.082851075618374 }
        };
    }

    //------------------------------------------------------------------------------
}    // namespace moris::mtk

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_12_HPP_ */
