/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_MTK_Integration_Coeffs_Bar_6.hpp
 *
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_32_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_32_HPP_

// MRS/COR/src
#include "moris_typedefs.hpp"
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
                Integration_Order::BAR_32 >::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_32 >::get_number_of_points()
        {
            return 32;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_32 >::get_points( Matrix< DDRMat > &aIntegrationPoints )
        {
            aIntegrationPoints = {
                { -0.9972638618494816,
                        -0.9856115115452684,
                        -0.9647622555875064,
                        -0.9349060759377397,
                        -0.8963211557660521,
                        -0.8493676137325700,
                        -0.7944837959679424,
                        -0.7321821187402897,
                        -0.6630442669302152,
                        -0.5877157572407623,
                        -0.5068999089322294,
                        -0.4213512761306353,
                        -0.3318686022821277,
                        -0.2392873622521371,
                        -0.1444719615827965,
                        -0.0483076656877383,
                        0.0483076656877383,
                        0.1444719615827965,
                        0.2392873622521371,
                        0.3318686022821277,
                        0.4213512761306353,
                        0.5068999089322294,
                        0.5877157572407623,
                        0.6630442669302152,
                        0.7321821187402897,
                        0.7944837959679424,
                        0.8493676137325700,
                        0.8963211557660521,
                        0.9349060759377397,
                        0.9647622555875064,
                        0.9856115115452684,
                        0.9972638618494816 }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_32 >::get_weights( Matrix< DDRMat > &aIntegrationWeights )
        {
            aIntegrationWeights = {
                {
                        0.0070186100094701,
                        0.0162743947309057,
                        0.0253920653092621,
                        0.0342738629130214,
                        0.0428358980222267,
                        0.0509980592623762,
                        0.0586840934785355,
                        0.0658222227763618,
                        0.0723457941088485,
                        0.0781938957870703,
                        0.0833119242269467,
                        0.0876520930044038,
                        0.0911738786957639,
                        0.0938443990808046,
                        0.0956387200792749,
                        0.0965400885147278,
                        0.0965400885147278,
                        0.0956387200792749,
                        0.0938443990808046,
                        0.0911738786957639,
                        0.0876520930044038,
                        0.0833119242269467,
                        0.0781938957870703,
                        0.0723457941088485,
                        0.0658222227763618,
                        0.0586840934785355,
                        0.0509980592623762,
                        0.0428358980222267,
                        0.0342738629130214,
                        0.0253920653092621,
                        0.0162743947309057,
                        0.0070186100094701,
                }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */
#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_BAR_32_HPP_ */
