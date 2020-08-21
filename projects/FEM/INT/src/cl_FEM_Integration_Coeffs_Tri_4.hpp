/*
 * cl_FEM_Integration_Coeffs_Tri_4.hpp
 *
 *  Created on: Aug 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_4_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_4_HPP_

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
        Integration_Order::TRI_4>::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_4>::get_number_of_points()
        {
            return 4;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_4>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                    {
                            0.333333333333333,

                            0.600000000000000,
                            0.200000000000000,
                            0.200000000000000
                    },
                    {
                            0.333333333333333,

                            0.200000000000000,
                            0.600000000000000,
                            0.200000000000000
                    },
                    {
                            0.333333333333333,

                            0.200000000000000,
                            0.200000000000000,
                            0.600000000000000
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_4 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            -0.562500000000000,

                            0.520833333333333,
                            0.520833333333333,
                            0.520833333333333
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_4_HPP_ */
