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

            // number of points in space
            uint mNumOfSpacePoints;

            // number of points in time
            uint mNumOfTimePoints;

            // matrix with space points
            Matrix< DDRMat > mSpacePoints;

            // matrix with time points
            Matrix< DDRMat > mTimePoints;

            // matrix with space weights
            Matrix< DDRMat > mSpaceWeights;

            // matrix with time weights
            Matrix< DDRMat > mTimeWeights;

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
            /**
             * get the number of integration points
             **/
            uint get_number_of_points();

//------------------------------------------------------------------------------
            /**
             * get the integration points
             **/
            void get_points( Matrix< DDRMat > & aIntegrationPoints );

//------------------------------------------------------------------------------
            /**
             * get the integration point weights
             **/
            void get_weights( Matrix< DDRMat > & aIntegrationWeights );

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATOR_HPP_ */
