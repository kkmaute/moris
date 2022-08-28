/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_Diffusion_Virtual_Work_Ghost_HPP_
#define SRC_FEM_CL_FEM_IWG_Diffusion_Virtual_Work_Ghost_HPP_

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

        class IWG_Diffusion_Virtual_Work_Ghost : public IWG
        {
                //------------------------------------------------------------------------------
            public:

                enum class IWG_Constitutive_Type
                {
                    DIFF_LIN_ISO,
                    MAX_ENUM
                };

                enum class IWG_Stabilization_Type
                {
                        GHOST_VW,
                        MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Diffusion_Virtual_Work_Ghost();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Diffusion_Virtual_Work_Ghost(){};

                //------------------------------------------------------------------------------
                /**
                 * compute the residual
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_jacobian( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the residual and the jacobian
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_jacobian_and_residual( real aWStar );

                //------------------------------------------------------------------------------
                /**
                 * compute the derivative of the residual wrt design variables
                 * @param[ in ] aWStar weight associated to the evaluation point
                 */
                void compute_dRdp( real aWStar );

            private:
                //------------------------------------------------------------------------------
                /**
                 * method to assemble "normal matrix" from normal vector needed for
                 * 2nd and 3rd order Ghost formulations
                 * @param[ in ] aFlatNormal flattened normal
                 * @param[ in ] aOrder      order of derivatives and ghost formulation
                 */
                void get_flat_normal_matrix(
                                Matrix< DDRMat > & aFlatNormal,
                                uint               aOrder );

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_Diffusion_Virtual_Work_Ghost_HPP_ */

