/*
 * cl_MTK_Integration_Coeffs_Tri_7.hpp
 *
 *  Created on: Jul 14, 2020
 *      Author: noel
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_7_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_7_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//MTK/src
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
        Integration_Order::TRI_7>::get_number_of_dimensions()
        {
            return 3;
        }
        //------------------------------------------------------------------------------


        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_7>::get_number_of_points()
        {
            return 7;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_7>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints = {
                    {
                            0.333333333333333,
                            0.797426985353087, 0.101286507323456, 0.101286507323456,
                            0.059715871789770, 0.470142064105115, 0.470142064105115
                    },
                    {
                            0.333333333333333,
                            0.101286507323456, 0.797426985353087, 0.101286507323456,
                            0.470142064105115, 0.059715871789770, 0.470142064105115
                    },
                    {
                            0.333333333333333,
                            0.101286507323456, 0.101286507323456, 0.797426985353087,
                            0.470142064105115, 0.470142064105115, 0.059715871789770
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_7>::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            0.225000000000000,
                            0.125939180544827, 0.125939180544827, 0.125939180544827,
                            0.132394152788506, 0.132394152788506, 0.132394152788506
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */

#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_TRI_7_HPP_ */