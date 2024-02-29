/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Interface_SLM_Mixed.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_INTERFACE_SLM_MIXED_HPP_
#define SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_INTERFACE_SLM_MIXED_HPP_

#include <map>

#include "moris_typedefs.hpp"                     //MRS/COR/src
#include "cl_Vector.hpp"                          //MRS/CNT/src

#include "cl_Matrix.hpp"                    //LINALG/src
#include "linalg_typedefs.hpp"              //LINALG/src

#include "cl_FEM_Field_Interpolator.hpp"    //FEM/INT/src
#include "cl_FEM_IWG.hpp"                   //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Isotropic_Struc_Linear_Interface_SLM_Mixed : public IWG
        {

            public:

                enum class IWG_Property_Type
                {
                	YOUNGS_MODULUS_LEADER,
					YOUNGS_MODULUS_FOLLOWER,
                    MAX_ENUM
                };

                enum class IWG_Constitutive_Type
                {
                        MAX_ENUM
                };

                enum class IWG_Stabilization_Type
                {
                        MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * @param[ in ] aBeta +1 or -1 for symmetric/unsymmetric symmetric Nitsche
                 */
                IWG_Isotropic_Struc_Linear_Interface_SLM_Mixed();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Isotropic_Struc_Linear_Interface_SLM_Mixed(){};

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

                //------------------------------------------------------------------------------
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_ISOTROPIC_STRUC_LINEAR_INTERFACE_HPP_ */

