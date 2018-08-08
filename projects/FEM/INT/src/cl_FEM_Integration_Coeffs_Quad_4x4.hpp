/*
 * cl_FEM_Integration_Coeffs_Quad_4x4.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_4X4_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_4X4_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src
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
            Integration_Order::QUAD_4x4>::get_number_of_dimensions()
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_4x4>::get_number_of_points()
            {
                return 16;
            }

//------------------------------------------------------------------------------

        template<>
        Mat< real >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_4x4>::get_points()
        {
            Mat< real > aIntegrationPoints =
            {
                {  -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                    0.339981043584856,
                    0.861136311594053
                 },
                 { -0.861136311594053,
                   -0.861136311594053,
                   -0.861136311594053,
                   -0.861136311594053,
                   -0.339981043584856,
                   -0.339981043584856,
                   -0.339981043584856,
                   -0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.339981043584856,
                    0.861136311594053,
                    0.861136311594053,
                    0.861136311594053,
                    0.861136311594053
                   }
            };

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Mat< real >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_4x4 >::get_weights()
            {
                Mat< real > aWeights =
                {
                   { 0.121002993285602,
                     0.226851851851852,
                     0.226851851851852,
                     0.121002993285602,
                     0.226851851851852,
                     0.425293303010694,
                     0.425293303010694,
                     0.226851851851852,
                     0.226851851851852,
                     0.425293303010694,
                     0.425293303010694,
                     0.226851851851852,
                     0.121002993285602,
                     0.226851851851852,
                     0.226851851851852,
                     0.121002993285602
                   }
                };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_4X4_HPP_ */
