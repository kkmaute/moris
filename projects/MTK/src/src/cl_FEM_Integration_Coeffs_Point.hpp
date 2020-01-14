/*
 * cl_FEM_Integration_Coeffs_Point.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_

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
        Integration_Order::POINT>::get_number_of_dimensions()
        {
            return 0;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::POINT>::get_number_of_points()
        {
                return 1;
        }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs< Integration_Type::GAUSS,
                            Integration_Order::POINT >::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints.set_size( 1, 1, 0.0 );
        }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::POINT >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights.set_size( 1, 1, 1.0 );
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_POINT_HPP_ */
