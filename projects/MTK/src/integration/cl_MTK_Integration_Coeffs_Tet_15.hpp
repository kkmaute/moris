/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Tet_15.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_15_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_15_HPP_

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
                Integration_Order::TET_15 >::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_15 >::get_number_of_points()
        {
            return 15;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_15 >::get_points( Matrix< DDRMat >& aIntegrationPoints )
        {
            aIntegrationPoints = {
                { //
                        0.2500000000000000,
                        0.0000000000000000,
                        0.3333333333333333,
                        0.3333333333333333,
                        0.3333333333333333,
                        0.7272727272727273,
                        0.0909090909090909,
                        0.0909090909090909,
                        0.0909090909090909,
                        0.4334498464263357,
                        0.0665501535736643,
                        0.0665501535736643,
                        0.0665501535736643,
                        0.4334498464263357,
                        0.4334498464263357 },
                { //
                        0.2500000000000000,
                        0.3333333333333333,
                        0.3333333333333333,
                        0.3333333333333333,
                        0.0000000000000000,
                        0.0909090909090909,
                        0.0909090909090909,
                        0.0909090909090909,
                        0.7272727272727273,
                        0.0665501535736643,
                        0.4334498464263357,
                        0.0665501535736643,
                        0.4334498464263357,
                        0.0665501535736643,
                        0.4334498464263357 },
                { //
                        0.2500000000000000,
                        0.3333333333333333,
                        0.0000000000000000,
                        0.3333333333333333,
                        0.3333333333333333,
                        0.0909090909090910,
                        0.7272727272727273,
                        0.0909090909090913,
                        0.0909090909090910,
                        0.4334498464263360,
                        0.4334498464263360,
                        0.4334498464263360,
                        0.0665501535736640,
                        0.0665501535736640,
                        0.0665501535736640 },
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TET_15 >::get_weights( Matrix< DDRMat >& aIntegrationWeights )
        {
            aIntegrationWeights = {
                { //
                        0.1817020685825351,
                        0.0361607142857143,
                        0.0361607142857143,
                        0.0361607142857143,
                        0.0361607142857143,
                        0.0698714945161738,
                        0.0698714945161738,
                        0.0698714945161738,
                        0.0698714945161738,
                        0.0656948493683187,
                        0.0656948493683187,
                        0.0656948493683187,
                        0.0656948493683187,
                        0.0656948493683187,
                        0.0656948493683187 }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TET_15_HPP_ */
