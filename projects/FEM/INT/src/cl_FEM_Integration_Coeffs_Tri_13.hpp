/*
 * cl_FEM_Integration_Coeffs_Tri_13.hpp
 *
 *  Created on: Aug 20, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_13_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_13_HPP_

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
        Integration_Order::TRI_13>::get_number_of_dimensions()
        {
            return 3;
        }

        //------------------------------------------------------------------------------

        template<>
        uint
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_13>::get_number_of_points()
        {
            return 13;
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_13>::get_points( Matrix< DDRMat > & aIntegrationPoints )
        {
            aIntegrationPoints =
            {
                    {
                            0.333333333333333,

                            0.479308067841920,
                            0.260345966079040,
                            0.260345966079040,

                            0.869739794195568,
                            0.065130102902216,
                            0.065130102902216,

                            0.048690315425316,
                            0.048690315425316,
                            0.312865496004874,
                            0.312865496004874,
                            0.638444188569810,
                            0.638444188569810
                    },
                    {
                            0.333333333333333,

                            0.260345966079040,
                            0.479308067841920,
                            0.260345966079040,

                            0.065130102902216,
                            0.869739794195568,
                            0.065130102902216,

                            0.312865496004874,
                            0.638444188569810,
                            0.048690315425316,
                            0.638444188569810,
                            0.048690315425316,
                            0.312865496004874
                    },
                    {
                            0.333333333333333,

                            0.260345966079040,
                            0.260345966079040,
                            0.479308067841920,

                            0.065130102902216,
                            0.065130102902216,
                            0.869739794195568,

                            0.638444188569810,
                            0.312865496004874,
                            0.638444188569810,
                            0.048690315425316,
                            0.312865496004874,
                            0.048690315425316
                    }
            };
        }

        //------------------------------------------------------------------------------

        template<>
        void
        Integration_Coeffs<
        Integration_Type::GAUSS,
        Integration_Order::TRI_13 >::get_weights( Matrix< DDRMat > & aIntegrationWeights )
        {
            aIntegrationWeights =
            {
                    {
                            -0.149570044467682,

                            0.175615257433208,
                            0.175615257433208,
                            0.175615257433208,

                            0.053347235608838,
                            0.053347235608838,
                            0.053347235608838,

                            0.077113760890257,
                            0.077113760890257,
                            0.077113760890257,
                            0.077113760890257,
                            0.077113760890257,
                            0.077113760890257
                    }
            };
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_TRI_13_HPP_ */
