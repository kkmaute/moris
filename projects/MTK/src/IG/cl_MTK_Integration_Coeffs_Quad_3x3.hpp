/*
 * cl_MTK_Integration_Coeffs_Quad_3x3.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_3X3_HPP_
#define SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_3X3_HPP_

#include "cl_MTK_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
#include "cl_MTK_Enums.hpp" //MTK/src

namespace moris
{
    namespace mtk
    {

//------------------------------------------------------------------------------

    template<>
    uint
    Integration_Coeffs<
    Integration_Type::GAUSS,
    Integration_Order::QUAD_3x3>::get_number_of_dimensions()
    {
        return 2;
    }


//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_3x3>::get_number_of_points()
            {
                return 9;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_3x3>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints = {
                    { -0.774596669241483,
                       0.774596669241483,
                       0.774596669241483,
                      -0.774596669241483,
                       0.000000000000000,
                       0.774596669241483,
                       0.000000000000000,
                      -0.774596669241483,
                       0.000000000000000 },
                    { -0.774596669241483,
                      -0.774596669241483,
                       0.774596669241483,
                       0.774596669241483,
                      -0.774596669241483,
                       0.000000000000000,
                       0.774596669241483,
                       0.000000000000000,
                       0.000000000000000 } };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_3x3 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights = {
                      { 0.308641975308642,
                        0.308641975308642,
                        0.308641975308642,
                        0.308641975308642,
                        0.493827160493827,
                        0.493827160493827,
                        0.493827160493827,
                        0.493827160493827,
                        0.790123456790123
                      } };
            }

//------------------------------------------------------------------------------
    } /* namespace mtk */
} /* namespace moris */


#endif /* SRC_MTK_CL_MTK_INTEGRATION_COEFFS_QUAD_3X3_HPP_ */