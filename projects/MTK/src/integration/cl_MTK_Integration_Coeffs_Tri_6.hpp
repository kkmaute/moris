/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tri_6.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_6_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_6_HPP_

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
                Integration_Order::TRI_6 >::get_number_of_dimensions()
        {
            return 2;
        }
        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_6 >::get_number_of_points()
        {
            return 6;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_6 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { 0.816847572980458,
                        0.091576213509771,
                        0.091576213509771,
                        0.108103018168070,
                        0.445948490915965,
                        0.445948490915965 },
                { 0.091576213509771,
                        0.816847572980458,
                        0.091576213509771,
                        0.445948490915965,
                        0.108103018168070,
                        0.445948490915965 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_6 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { 0.109951743655322,
                        0.109951743655322,
                        0.109951743655322,
                        0.223381589678011,
                        0.223381589678011,
                        0.223381589678011 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_6_HPP_ */

