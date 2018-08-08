/*
 * cl_FEM_Integrator.cpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BASE_HPP_
#define SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BASE_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        class Integration_Coeffs_Base
        {

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            /* trivial constructor */
            Integration_Coeffs_Base(){};

//------------------------------------------------------------------------------

            /* trivial destructor */
            virtual
            ~Integration_Coeffs_Base(){};

//------------------------------------------------------------------------------

            /**
             * returns the number of dimensions
             */
            virtual uint
            get_number_of_dimensions() = 0;

//------------------------------------------------------------------------------

            /**
             * returns the number of points
             */
            virtual uint
            get_number_of_points() = 0;

//------------------------------------------------------------------------------

            /**
             * returns the integration weights
             *
             * @param[ out ] aIntegrationWeights
             */
            virtual Mat< real >
            get_weights() = 0;

//------------------------------------------------------------------------------

            /**
             * writes the integration points into given Mat
             *
             * @param[ out ] aIntegrationPoints
             */
            virtual Mat< real >
            get_points() = 0;

//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATION_COEFFS_BASE_HPP_ */
