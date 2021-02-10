/*
 * cl_FEM_Integration_Coeffs_Bar_1.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_1_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_1_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_FEM_Enums.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::BAR_1>::get_number_of_dimensions()
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_1>::get_number_of_points()
        {
                return 1;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs< Integration_Type::GAUSS,
                            Integration_Order::BAR_1 >::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_1 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights.set_size( 1, 1, 2.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_1_HPP_ */
