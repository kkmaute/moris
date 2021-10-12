/*
 * cl_FEM_IQI_Factory.cpp
 *
 *  Created on: Dec 5, 2019
 *      Author: noel
 */

#include <memory>
#include "assert.hpp"
//FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IQI_Volume.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"
#include "cl_FEM_IQI_Stress.hpp"
#include "cl_FEM_IQI_Dof.hpp"
#include "cl_FEM_IQI_Max_Dof.hpp"
#include "cl_FEM_IQI_Property.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"
#include "cl_FEM_IQI_H1_Error_Analytic.hpp"
#include "cl_FEM_IQI_J_Integral.hpp"
#include "cl_FEM_IQI_Volume_Fraction.hpp"
#include "cl_FEM_IQI_2D_Drag_Lift_Coefficient.hpp"
#include "cl_FEM_IQI_H1_Error.hpp"
#include "cl_FEM_IQI_Latent_Heat_Absorption.hpp"
#include "cl_FEM_IQI_Max_Stress.hpp"
#include "cl_FEM_IQI_Turbulent_Kinematic_Viscosity.hpp"
#include "cl_FEM_IQI_Total_Pressure.hpp"
#include "cl_FEM_IQI_Mass_Flow.hpp"
#include "cl_FEM_IQI_Stabilization.hpp"
#include "cl_FEM_IQI_Homogenized_Constitutive.hpp"
#include "cl_FEM_IQI_Heat_Method_Penalty.hpp"
#include "cl_FEM_IQI_Thermal_Energy_Convective_Flux.hpp"
#include "cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp"
#include "cl_FEM_IQI_Advection_Strong_Residual.hpp"
#include "cl_FEM_IQI_Zienkiewicz_Zhu.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        std::shared_ptr< IQI > IQI_Factory::create_IQI( IQI_Type aIQIType )
        {
            switch( aIQIType )
            {
                case IQI_Type::DOF :
                    return std::make_shared< IQI_Dof >();

                case IQI_Type::MAX_DOF :
                    return std::make_shared< IQI_Max_Dof >();

                case IQI_Type::PROPERTY :
                    return std::make_shared< IQI_Property >();

                case IQI_Type::STABILIZATION :
                	return std::make_shared< IQI_Stabilization >();

                case IQI_Type::VOLUME :
                    return std::make_shared< IQI_Volume >();

                case IQI_Type::VOLUME_FRACTION :
                    return std::make_shared< IQI_Volume_Fraction >();

                case IQI_Type::STRAIN_ENERGY :
                    return std::make_shared< IQI_Strain_Energy >();

                case IQI_Type::NORMAL_STRESS :
                    return std::make_shared< IQI_Stress >( Stress_Type::NORMAL_STRESS );
                case IQI_Type::SHEAR_STRESS :
                    return std::make_shared< IQI_Stress >( Stress_Type::SHEAR_STRESS );
                case IQI_Type::VON_MISES_STRESS :
                    return std::make_shared< IQI_Stress >( Stress_Type::VON_MISES_STRESS );
                case IQI_Type::PRINCIPAL_STRESS :
                    return std::make_shared< IQI_Stress >( Stress_Type::PRINCIPAL_STRESS );

                case IQI_Type::MAX_NORMAL_STRESS :
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::NORMAL_STRESS );
                case IQI_Type::MAX_SHEAR_STRESS :
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::SHEAR_STRESS );
                case IQI_Type::MAX_VON_MISES_STRESS :
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::VON_MISES_STRESS );
                case IQI_Type::MAX_PRINCIPAL_STRESS :
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::PRINCIPAL_STRESS );

                case IQI_Type::L2_ERROR_ANALYTIC :
                    return std::make_shared< IQI_L2_Error_Analytic >();

                case IQI_Type::H1_ERROR_ANALYTIC :
                    return std::make_shared< IQI_H1_Error_Analytic >();

                case IQI_Type::H1_ERROR :
                    return std::make_shared< IQI_H1_Error >();

                case IQI_Type::J_INTEGRAL :
                    return std::make_shared< IQI_J_Integral >();

                case IQI_Type::DRAG_COEFF :
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( 1 );

                case IQI_Type::LIFT_COEFF :
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( -1 );

                case IQI_Type::TURBULENT_KINEMATIC_VISCOSITY :
                    return std::make_shared< IQI_Turbulent_Kinematic_Viscosity >();

                case IQI_Type::LATENT_HEAT_ABSORPTION :
                    return std::make_shared< IQI_Latent_Heat_Absorption >();

                case IQI_Type::TOTAL_PRESSURE :
                    return std::make_shared< IQI_Total_Pressure >();

                case IQI_Type::MASS_FLOW :
                    return std::make_shared< IQI_Mass_Flow >();

                case IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX :
                    return std::make_shared< IQI_Thermal_Energy_Convective_Flux >();

                case IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX :
                    return std::make_shared< IQI_Thermal_Energy_Diffusive_Flux >();

                case IQI_Type::HOMOGENIZED_CONSTITUTIVE:
                    return std::make_shared< IQI_Homogenized_Constitutive >();

                case IQI_Type::HEAT_METHOD_PENALTY:
                    return std::make_shared< IQI_Heat_Method_Penalty >();

                case IQI_Type::ADVECTION_STRONG_RESIDUAL:
                    return std::make_shared< IQI_Advection_Strong_Residual >();

                case IQI_Type::ZIENKIEWICZ_ZHU_VON_MISES_STRESS:
                    return std::make_shared< IQI_Zienkiewicz_Zhu >( Stress_Type::VON_MISES_STRESS );

                default:
                    MORIS_ERROR( false, " IQI_Factory::create_IQI - No IQI type specified. " );
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

