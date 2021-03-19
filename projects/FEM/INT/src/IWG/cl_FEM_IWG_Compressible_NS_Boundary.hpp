/*
 * cl_FEM_IWG_Compressible_NS_Boundary.hpp
 *
 *  Created on: Mar 16, 2021
 *      Author: wunsch
 */

#ifndef SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_
#define SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_

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

        class IWG_Compressible_NS_Boundary : public IWG
        {
                //------------------------------------------------------------------------------

            private:

                // default dof types (FIXME: only primitive vars for now)
                MSI::Dof_Type mDofDensity     = MSI::Dof_Type::RHO;
                MSI::Dof_Type mDofPressure    = MSI::Dof_Type::P;
                MSI::Dof_Type mDofVelocity    = MSI::Dof_Type::VX;
                MSI::Dof_Type mDofTemperature = MSI::Dof_Type::TEMP;

                // flag indicating what variable set is being used
                fem::Variable_Set mVariableSet = fem::Variable_Set::UNDEFINED;

                //------------------------------------------------------------------------------

            public:

                enum class IWG_Property_Type
                {
                    HEAT_FLUX,
                    TRACTION,
                    PRESSURE,
                    MAX_ENUM
                };

                // local material enums
                enum class IWG_Material_Type
                {
                    FLUID_MM,
                    MAX_ENUM
                };

                // local constitutive enums
                enum class IWG_Constitutive_Type
                {
                    FLUID_CM,
                    MAX_ENUM
                };

                //------------------------------------------------------------------------------
                /*
                 * constructor
                 */
                IWG_Compressible_NS_Boundary();

                //------------------------------------------------------------------------------
                /**
                 * trivial destructor
                 */
                ~IWG_Compressible_NS_Boundary(){};

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
                /**
                 * check that the set of residual DoF types is valid
                 * and supported by the implementation
                 */
                bool check_residual_dof_types();

                //------------------------------------------------------------------------------

        };
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_IWG_COMPRESSIBLE_NS_BOUNDARY_HPP_ */
