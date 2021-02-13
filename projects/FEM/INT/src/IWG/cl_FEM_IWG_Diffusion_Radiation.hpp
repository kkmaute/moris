/*
 * cl_FEM_IWG_Diffusion_Radiation.hpp
 *
 *  Created on: July 17, 2020
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_Diffusion_Radiation_HPP_
#define SRC_FEM_CL_FEM_IWG_Diffusion_Radiation_HPP_

#include <map>

#include "typedefs.hpp"                     //MRS/COR/src
#include "cl_Cell.hpp"                      //MRS/CON/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Diffusion_Radiation : public IWG
        {
            private:

                // Stefan-Bolzmann constant for black body radiation
                const real mStefanBoltzmannConst = 5.670374419e-08;

                //------------------------------------------------------------------------------
            public:
                enum class IWG_Property_Type
                {
                    EMISSIVITY,
                    AMBIENT_TEMP,
                    ABSOLUTE_ZERO,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IWG_Diffusion_Radiation();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Diffusion_Radiation(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_jacobian( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aWStar weight associated with evaluation point
                 */
                void compute_jacobian_and_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the residual wrt design variables
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dRdp( real aWStar );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_Diffusion_Radiation_HPP_ */