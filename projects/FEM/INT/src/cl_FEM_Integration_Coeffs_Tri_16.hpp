/*
 * cl_FEM_Integration_Coeffs_Tri_16.hpp
 *
 *  Created on: Aug 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_16_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_16_HPP_

//MRS/COR/src
#include "typedefs.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Integration_Coeffs.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_16>::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_16>::get_number_of_points()
        {
            return 16;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_16>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                    {
                            0.333333333333333,

                            0.081414823414554,
                            0.459292588292723,
                            0.459292588292723,

                            0.658861384496480,
                            0.170569307751760,
                            0.170569307751760,

                            0.898905543365938,
                            0.050547228317031,
                            0.050547228317031,

                            0.008394777409958,
                            0.008394777409958,
                            0.263112829634638,
                            0.263112829634638,
                            0.728492392955404,
                            0.728492392955404
                    },
                    {
                            0.333333333333333,

                            0.459292588292723,
                            0.081414823414554,
                            0.459292588292723,

                            0.170569307751760,
                            0.658861384496480,
                            0.170569307751760,

                            0.050547228317031,
                            0.898905543365938,
                            0.050547228317031,

                            0.263112829634638,
                            0.728492392955404,
                            0.008394777409958,
                            0.728492392955404,
                            0.008394777409958,
                            0.263112829634638
                    },
                    {
                            0.333333333333333,

                            0.459292588292723,
                            0.459292588292723,
                            0.081414823414554,

                            0.170569307751760,
                            0.170569307751760,
                            0.658861384496480,

                            0.050547228317031,
                            0.050547228317031,
                            0.898905543365938,

                            0.728492392955404,
                            0.263112829634638,
                            0.728492392955404,
                            0.008394777409958,
                            0.263112829634638,
                            0.008394777409958
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_16 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            0.144315607677787,

                            0.095091634267285,
                            0.095091634267285,
                            0.095091634267285,

                            0.103217370534718,
                            0.103217370534718,
                            0.103217370534718,

                            0.032458497623198,
                            0.032458497623198,
                            0.032458497623198,

                            0.027230314174435,
                            0.027230314174435,
                            0.027230314174435,
                            0.027230314174435,
                            0.027230314174435,
                            0.027230314174435
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_16_HPP_ */
