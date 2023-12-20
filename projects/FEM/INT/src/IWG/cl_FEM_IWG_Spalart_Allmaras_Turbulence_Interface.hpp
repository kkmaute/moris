/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Spalart_Allmaras_Turbulence_Interface.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_
#define SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_
//MRS/COR/src
#include <map>
#include "moris_typedefs.hpp"
#include "cl_Cell.hpp"
//LINALG/src
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
//FEM/INT/src
#include "cl_FEM_IWG.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        class IWG_Spalart_Allmaras_Turbulence_Interface : public IWG
        {

                //------------------------------------------------------------------------------
            private:

                // sint for symmetric/unsymmetric Nitsche formulation
                sint mBeta = 1.0;

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        SPALART_ALLMARAS_TURBULENCE,
                        MAX_ENUM
                };

                // local stabilization parameter enums
                enum class IWG_Stabilization_Type
                {
                        NITSCHE_INTERFACE,
                        MAX_ENUM
                };

            public:
                //------------------------------------------------------------------------------
                /*
                 * constructor
                 * @param[ in ] aBeta sint for symmetric/unsymmetric Nitsche formulation
                 */
                IWG_Spalart_Allmaras_Turbulence_Interface( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Spalart_Allmaras_Turbulence_Interface(){};

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
        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_SPALART_ALLMARAS_TURBULENCE_INTERFACE_HPP_ */

