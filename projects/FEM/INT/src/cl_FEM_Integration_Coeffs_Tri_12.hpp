/*
 * cl_FEM_Integration_Coeffs_Tri_12.hpp
 *
 *  Created on: Jul 15, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_12_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_12_HPP_

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
        Integration_Order::TRI_12>::get_number_of_dimensions()
        {
            return 3;
        }
        //------------------------------------------------------------------------------


        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_12>::get_number_of_points()
        {
            return 12;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_12>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints = {
                    {
                            0.873821971016996, 0.063089014491502, 0.063089014491502,
                            0.501426509658179, 0.249286745170910, 0.249286745170910,
                            0.636502499121399, 0.636502499121399, 0.310352451033785, 0.310352451033785, 0.053145049844816, 0.053145049844816
                    },
                    {
                            0.063089014491502, 0.873821971016996, 0.063089014491502,
                            0.249286745170910, 0.501426509658179, 0.249286745170910,
                            0.310352451033785, 0.053145049844816, 0.636502499121399, 0.053145049844816, 0.310352451033785, 0.636502499121399
                    },
                    {
                            0.063089014491502, 0.063089014491502, 0.873821971016996,
                            0.249286745170910, 0.249286745170910, 0.501426509658179,
                            0.053145049844816, 0.310352451033785, 0.053145049844816, 0.636502499121399, 0.636502499121399, 0.310352451033785
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_12>::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            0.050844906370207, 0.050844906370207, 0.050844906370207,
                            0.116786275726379, 0.116786275726379, 0.116786275726379,
                            0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374, 0.082851075618374
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_12_HPP_ */
