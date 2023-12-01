/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tet_11.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_11_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_11_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "typedefs.hpp"           //MRS/COR/src
#include "cl_Matrix.hpp"          //LNA/src
#include "linalg_typedefs.hpp"    //LNA/src
#include "cl_MTK_Enums.hpp"       //MTK/src
#include "op_times.hpp"           //MTK/src

namespace moris
{
    namespace mtk
    {

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11 >::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11 >::get_number_of_points()
        {
            return 11;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { //
                        0.7857142857142860,
                        0.0714285714285714,
                        0.0714285714285714,
                        0.0714285714285714,
                        0.3994035761667990,
                        0.3994035761667990,
                        0.3994035761667990,
                        0.1005964238332010,
                        0.1005964238332010,
                        0.1005964238332010,
                        0.2500000000000000 },
                { //
                        0.0714285714285714,
                        0.7857142857142860,
                        0.0714285714285714,
                        0.0714285714285714,
                        0.3994035761667990,
                        0.1005964238332010,
                        0.1005964238332010,
                        0.3994035761667990,
                        0.3994035761667990,
                        0.1005964238332010,
                        0.2500000000000000 },
                { //
                        0.0714285714285714,
                        0.0714285714285714,
                        0.0714285714285714,
                        0.7857142857142860,
                        0.1005964238332010,
                        0.1005964238332010,
                        0.3994035761667990,
                        0.1005964238332010,
                        0.3994035761667990,
                        0.3994035761667990,
                        0.2500000000000000 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_11 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { //
                        0.0076222222222222,
                        0.0076222222222222,
                        0.0076222222222222,
                        0.0076222222222222,
                        0.0248888888888889,
                        0.0248888888888889,
                        0.0248888888888889,
                        0.0248888888888889,
                        0.0248888888888889,
                        0.0248888888888889,
                        -0.0131555555555556 }
            };

            aIntegrationWeights = 6.0 * aIntegrationWeights;
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_11_HPP_ */
