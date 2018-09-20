/*
 * cl_FEM_Integration_Coeffs_Quad_2x2.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_2X2_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_2X2_HPP_

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::QUAD_2x2>::get_number_of_dimensions()
        {
            return 2;
        }

//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::QUAD_2x2>::get_number_of_points()
            {
                return 4;
            }

//------------------------------------------------------------------------------

        template<>
        Matrix< DDRMat >
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_2x2>::get_points()
        {
            Matrix< DDRMat > aIntegrationPoints = {
                    { -0.577350269189626,
                       0.577350269189626,
                       0.577350269189626,
                      -0.577350269189626 },
                    { -0.577350269189626,
                      -0.577350269189626,
                       0.577350269189626,
                       0.577350269189626 } };

            return aIntegrationPoints;
          }

//------------------------------------------------------------------------------

            template<>
            Matrix< DDRMat >
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::QUAD_2x2 >::get_weights()
            {
                Matrix< DDRMat > aWeights = {
                        { 1.000000000000000,
                          1.000000000000000,
                          1.000000000000000,
                          1.000000000000000} };

                return aWeights;
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_QUAD_2X2_HPP_ */
