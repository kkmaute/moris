/*
 * cl_FEM_IWG_Compressible_NS_Velocity_Dirichlet_Nitsche.hpp
 *
 *  Created on: Dec 13, 2020
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_

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

        class IWG_Compressible_NS_Velocity_Dirichlet_Nitsche : public IWG
        {

                //------------------------------------------------------------------------------
            public:

                // sign for symmetric/unsymmetric Nitsche
                sint mBeta;

                // default dof types
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                enum class IWG_Property_Type
                {
                        PRESCRIBED_VALUE,
                        SELECT,
                        MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                        FLUID,
                        MAX_ENUM
                };

                // local stabilization enums
                enum class IWG_Stabilization_Type
                {
                        NITSCHE_PENALTY_PARAMETER,
                        MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 *  constructor
                 */
                IWG_Compressible_NS_Velocity_Dirichlet_Nitsche( sint aBeta );

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Velocity_Dirichlet_Nitsche(){};

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

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_VELOCITY_DIRICHLET_NITSCHE_HPP_ */