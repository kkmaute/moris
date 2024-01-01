/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_Diffusion_Dirichlet_Nitsche_HPP_
#define SRC_FEM_CL_FEM_IWG_Diffusion_Dirichlet_Nitsche_HPP_

#include <map>

#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"     //MRS/CNT/src

#include "cl_Matrix.hpp"          //LINALG/src
#include "linalg_typedefs.hpp"    //LINALG/src

#include "cl_FEM_IWG.hpp"    //FEM/INT/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Diffusion_Dirichlet_Nitsche : public IWG
        {

            //------------------------------------------------------------------------------

          public:
            // sign for symmetric/unsymmetric Nitsche
            sint mBeta;

            enum class IWG_Property_Type
            {
                DIRICHLET,
                SELECT,
                THICKNESS,
                MAX_ENUM
            };

            enum class IWG_Constitutive_Type
            {
                DIFF_LIN_ISO,
                MAX_ENUM
            };

            enum class IWG_Stabilization_Type
            {
                DIRICHLET_NITSCHE,
                MAX_ENUM
            };

            //------------------------------------------------------------------------------
            /*
             * constructor
             * @param[ in ] aBeta sint for symmetric/unsymmetric Nitsche
             */
            IWG_Diffusion_Dirichlet_Nitsche( sint aBeta );

            //------------------------------------------------------------------------------
            /**
             * trivial destructor
             */
            ~IWG_Diffusion_Dirichlet_Nitsche(){};

            //------------------------------------------------------------------------------
            /**
             * computes the residual
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_residual( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * computes the jacobian
             * @param[ in ] aWStar weight associated to the evaluation point
             */
            void compute_jacobian( real aWStar );

            //------------------------------------------------------------------------------
            /**
             * computes the residual and the jacobian
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

#endif /* SRC_FEM_CL_FEM_IWG_Diffusion_Dirichlet_Nitsche_HPP_ */
