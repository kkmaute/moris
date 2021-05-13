/*
 * cl_FEM_IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche.hpp
 *
 *  Created on: Mar 25, 2020
 *      Author: noel
 */

#ifndef SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_

#include <map>
//MRS/COR/src
#include "typedefs.hpp"
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

        class IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // sign for symmetric/unsymmetric Nitsche
                sint mBeta;

                enum class IWG_Property_Type
                {
                        DIRICHLET,
                        SELECT,
                        UPWIND,
                        MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        FLUID_INCOMPRESSIBLE,
                        MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                        VELOCITY_DIRICHLET_NITSCHE,
                        MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Incompressible_NS_Velocity_Dirichlet_Nitsche(){};

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

#endif /* SRC_FEM_CL_FEM_IWG_INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_ */
