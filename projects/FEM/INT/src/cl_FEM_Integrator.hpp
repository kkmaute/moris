/*
 * cl_FEM_Integrator.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATOR_HPP_
#define SRC_FEM_CL_FEM_INTEGRATOR_HPP_

#include "typedefs.hpp"                  //MRS/COR/src
#include "cl_Matrix.hpp"                 //LNA/src
#include "cl_FEM_Enums.hpp"              //FEM/INT/src
#include "cl_FEM_Integration_Rule.hpp"   //FEM/INT/src
#include "cl_FEM_Integration_Coeffs.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Integrator
        {
            // pointer to space rule, if specified
            Integration_Coeffs_Base * mSpaceCoeffs      = nullptr;

            // pointer to time rule, if specified
            Integration_Coeffs_Base * mTimeCoeffs       = nullptr;

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------
            /**
             * constructs an integrator from an integration rule
             **/
            Integrator( const Integration_Rule & aIntegrationRule );

//------------------------------------------------------------------------------
            /**
             * destructor
             **/
            ~Integrator();

//------------------------------------------------------------------------------
//            /**
//             * get the number of dimensions
//             **/
//            uint get_number_of_dimensions();

//------------------------------------------------------------------------------
            /**
             * get the number of integration points
             **/
            uint get_number_of_points();

//------------------------------------------------------------------------------
            /**
             * get the integration points
             **/
            Matrix< DDRMat > get_points();

//------------------------------------------------------------------------------
            /**
             * get the integration point weights
             **/
            Matrix< DDRMat > get_weights();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATOR_HPP_ */
