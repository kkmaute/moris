/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tet_5.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_5_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_5_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "typedefs.hpp"           //MRS/COR/src
#include "cl_Matrix.hpp"          //LNA/src
#include "linalg_typedefs.hpp"    //LNA/src
#include "cl_MTK_Enums.hpp"       //MTK/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------
        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5 >::get_number_of_dimensions()
        {
            return 3;
        }
        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5 >::get_number_of_points()
        {
            return 5;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { //
                        0.500000000000000,
                        0.166666666666667,
                        0.166666666666667,
                        0.166666666666667,
                        0.250000000000000 },
                { //
                        0.166666666666667,
                        0.500000000000000,
                        0.166666666666667,
                        0.166666666666667,
                        0.250000000000000 },
                { //
                        0.166666666666667,
                        0.166666666666667,
                        0.166666666666667,
                        0.500000000000000,
                        0.250000000000000 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_5 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { //
                        0.450000000000000,
                        0.450000000000000,
                        0.450000000000000,
                        0.450000000000000,
                        -0.800000000000000 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_5_HPP_ */
