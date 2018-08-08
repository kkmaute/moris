/*
 * cl_FEM_Integrator.hpp
 *
 *  Created on: Jul 19, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_INTEGRATOR_HPP_
#define SRC_FEM_CL_FEM_INTEGRATOR_HPP_

#include "typedefs.hpp" //MRS/COR/src
#include "cl_Mat.hpp" //LNA/src
#include "cl_FEM_Enums.hpp" //FEM/INT/src
#include "cl_FEM_Integration_Rule.hpp" //FEM/INT/src
#include "cl_FEM_Integration_Coeffs.hpp" //FEM/INT/src

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        class Integrator
        {
            //! pointer to space rule, if specified
            Integration_Coeffs_Base * mSpaceCoeffs      = nullptr;

            //! pointer to time rule, if specified
            Integration_Coeffs_Base * mTimeCoeffs       = nullptr;

            //! pointer to space time rule, if specified
            Integration_Coeffs_Base * mSpaceTimeCoeffs  = nullptr;

            //! pointer to dimensions
            uint
            ( Integrator:: * mGetNumberOfDimensions )();

            //! pointer to get_number_of_points function
            uint
            ( Integrator:: * mGetNumberOfPoints )();

            //! pointer to get_points function
            Mat< real >
            ( Integrator:: * mGetPoints )();

            //! pointer to get_weights_function
            Mat< real >
            ( Integrator:: * mGetWeights )();

//------------------------------------------------------------------------------
        public :
//------------------------------------------------------------------------------

            Integrator( const Integration_Rule & aIntegrationRule );

//------------------------------------------------------------------------------

            ~Integrator();

//------------------------------------------------------------------------------

            uint
            get_number_of_dimensions();

//------------------------------------------------------------------------------

            uint
            get_number_of_points();

//------------------------------------------------------------------------------

            Mat< real >
            get_points();

//------------------------------------------------------------------------------

            Mat< real >
            get_weights();

//------------------------------------------------------------------------------
        private :
//------------------------------------------------------------------------------

            uint
            get_number_of_dimensions_space_and_time();

//------------------------------------------------------------------------------

            uint
            get_number_of_points_space_and_time();

//------------------------------------------------------------------------------

            Mat< real >
            get_points_space_and_time();

//------------------------------------------------------------------------------

            Mat< real >
            get_weights_space_and_time();

//------------------------------------------------------------------------------

            uint
            get_number_of_dimensions_spacetime();

//------------------------------------------------------------------------------

            uint
            get_number_of_points_spacetime();

//------------------------------------------------------------------------------

            Mat< real >
            get_points_spacetime();

//------------------------------------------------------------------------------

            Mat< real >
            get_weights_spacetime();

//------------------------------------------------------------------------------
        };
//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_INTEGRATOR_HPP_ */
