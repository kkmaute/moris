/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tet_4.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_4_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_4_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "moris_typedefs.hpp"           //MRS/COR/src
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
                Integration_Order::TET_4 >::get_number_of_dimensions()
        {
            return 3;
        }
        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_4 >::get_number_of_points()
        {
            return 4;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_4 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { //
                        0.585410196624968,
                        0.138196601125010,
                        0.138196601125010,
                        0.138196601125010 },
                { //
                        0.138196601125010,
                        0.585410196624968,
                        0.138196601125010,
                        0.138196601125010 },
                { //
                        0.138196601125010,
                        0.138196601125010,
                        0.138196601125010,
                        0.585410196624968 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_4 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { //
                        0.250000000000000,
                        0.250000000000000,
                        0.250000000000000,
                        0.250000000000000 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_4_HPP_ */
