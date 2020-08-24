/*
 * cl_FEM_Integration_Coeffs_Tri_1.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_6_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_6_HPP_

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
            Integration_Order::TRI_6>::get_number_of_dimensions()
        {
            return 3;
        }
//------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
            Integration_Type::GAUSS,
            Integration_Order::TRI_6>::get_number_of_points()
            {
                return 6;
            }

//------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_6>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                {
                    0.816847572980458,
                    0.091576213509771,
                    0.091576213509771,
                    0.108103018168070,
                    0.445948490915965,
                    0.445948490915965
                },
                {
                    0.091576213509771,
                    0.816847572980458,
                    0.091576213509771,
                    0.445948490915965,
                    0.108103018168070,
                    0.445948490915965
                },
                {
                    0.091576213509771,
                    0.091576213509771,
                    0.816847572980458,
                    0.445948490915965,
                    0.445948490915965,
                    0.108103018168070
                }
            };
          }

//------------------------------------------------------------------------------

            template<>
            void
            Integration_Coeffs<
                Integration_Type::GAUSS,
                Integration_Order::TRI_6 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
            {
                aIntegrationWeights =
                {
                    {
                        0.109951743655322,
                        0.109951743655322,
                        0.109951743655322,
                        0.223381589678011,
                        0.223381589678011,
                        0.223381589678011
                     }
                };
            }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_6_HPP_ */
