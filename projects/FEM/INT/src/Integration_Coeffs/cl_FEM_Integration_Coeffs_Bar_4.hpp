/*
 * cl_FEM_Integration_Coeffs_Bar_2.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_4_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_4_HPP_

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
        Integration_Order::BAR_4>::get_number_of_dimensions()
        {
            return 1;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::BAR_4>::get_number_of_points()
            {
                return 4;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_4>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                 {
                     -8.611363115940526e-01,
                     -3.399810435848563e-01,
                      3.399810435848563e-01,
                      8.611363115940526e-01
                 }

            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::BAR_4 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                    {
                        3.478548451374538e-01,
                        6.521451548625461e-01,
                        6.521451548625461e-01,
                        3.478548451374538e-01
                    }

                };
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BAR_4_HPP_ */
