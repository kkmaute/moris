/*
 * cl_FEM_Integration_Coeffs_Quad_3x3.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_3X3_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_3X3_HPP_

#include "cl_FEM_Integration_Coeffs.hpp"
#include "typedefs.hpp" //MRS/COR/src
#include "cl_Matrix.hpp" //LNA/src
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
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_3x3>::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints = {
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


            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_3x3 >::get_weights()
            {
                Matrix< DDRMat > aWeights = {
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

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_3X3_HPP_ */
