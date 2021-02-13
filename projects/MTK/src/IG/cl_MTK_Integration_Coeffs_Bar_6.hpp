/*
 * cl_FEM_Integration_Coeffs_Bar_6.hpp
 *
 *  Created on: Aug 21, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_6_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_6_HPP_

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
        Integration_Order::BAR_6>::get_number_of_dimensions()
        {
            return 1;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::BAR_6>::get_number_of_points()
        {
            return 6;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::BAR_6>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                    {
                            -0.932469514203152,
                            -0.661209386466265,
                            -0.238619186083197,
                            0.238619186083197,
                            0.661209386466265,
                            0.932469514203152
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::BAR_6 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            0.171324492379170,
                            0.360761573048139,
                            0.467913934572691,
                            0.467913934572691,
                            0.360761573048139,
                            0.171324492379170
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_6_HPP_ */
