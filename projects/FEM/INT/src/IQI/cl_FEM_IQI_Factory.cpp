/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Factory.cpp
 *
 */

#include <memory>
#include "assert.hpp"
// FEM/INT/src
#include "cl_FEM_IQI_Factory.hpp"
#include "cl_FEM_IQI_Volume.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"
#include "cl_FEM_IQI_Stress.hpp"
#include "cl_FEM_IQI_Gap.hpp"
#include "cl_FEM_IQI_Dof.hpp"
#include "cl_FEM_IQI_Eigen_Vector.hpp"
#include "cl_FEM_IQI_Eigen_Value.hpp"
#include "cl_FEM_IQI_ALM_Dof.hpp"
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
#include "cl_FEM_IQI_Turbulent_Dynamic_Viscosity.hpp"
#include "cl_FEM_IQI_Effective_Dynamic_Viscosity.hpp"
#include "cl_FEM_IQI_Effective_Conductivity.hpp"
#include "cl_FEM_IQI_Spalart_Allmaras_Coefficient.hpp"
#include "cl_FEM_IQI_Power_Dissipation.hpp"
#include "cl_FEM_IQI_Power_Dissipation_Bulk.hpp"
#include "cl_FEM_IQI_Total_Pressure.hpp"
#include "cl_FEM_IQI_Mass_Flow.hpp"
#include "cl_FEM_IQI_Stabilization.hpp"
#include "cl_FEM_IQI_Homogenized_Constitutive.hpp"
#include "cl_FEM_IQI_Heat_Method_Penalty.hpp"
#include "cl_FEM_IQI_Thermal_Energy_Convective_Flux.hpp"
#include "cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.hpp"
#include "cl_FEM_IQI_Advection_Strong_Residual.hpp"
#include "cl_FEM_IQI_Strong_Residual_SA.hpp"
#include "cl_FEM_IQI_Strong_Residual_Incompressible_NS.hpp"
#include "cl_FEM_IQI_SP_Crosswind_Incompressible_NS.hpp"
#include "cl_FEM_IQI_Res_Crosswind_Incompressible_NS.hpp"
#include "cl_FEM_IQI_Res_SUPG_Incompressible_NS.hpp"
#include "cl_FEM_IQI_SP_Crosswind_SA.hpp"
#include "cl_FEM_IQI_Zienkiewicz_Zhu.hpp"
#include "cl_FEM_IQI_Normal_Vector.hpp"
#include "cl_FEM_IQI_Jump_Dof.hpp"
#include "cl_FEM_IQI_Jump_Traction.hpp"
#include "cl_FEM_IQI_Traction.hpp"
#include "cl_FEM_IQI_Linear_Elasticity_Damage.hpp"
#include "cl_FEM_IQI_Contact_Pressure.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        std::shared_ptr< IQI >
        IQI_Factory::create_IQI( IQI_Type aIQIType )
        {
            switch ( aIQIType )
            {
                case IQI_Type::DOF:
                    return std::make_shared< IQI_Dof >();

                case IQI_Type::EIGEN_VECTOR:
                    return std::make_shared< IQI_Eigen_Vector >();

                case IQI_Type::EIGEN_VALUE:
                    return std::make_shared< IQI_Eigen_Value >();

                case IQI_Type::ALM_DOF:
                    return std::make_shared< IQI_ALM_Dof >();

                case IQI_Type::MAX_DOF:
                    return std::make_shared< IQI_Max_Dof >();

                case IQI_Type::PROPERTY:
                    return std::make_shared< IQI_Property >();

                case IQI_Type::STABILIZATION:
                    return std::make_shared< IQI_Stabilization >();

                case IQI_Type::VOLUME:
                    return std::make_shared< IQI_Volume >();

                case IQI_Type::VOLUME_FRACTION:
                    return std::make_shared< IQI_Volume_Fraction >();

                case IQI_Type::STRAIN_ENERGY:
                    return std::make_shared< IQI_Strain_Energy >();

                case IQI_Type::NORMAL_STRESS:
                    return std::make_shared< IQI_Stress >( Stress_Type::NORMAL_STRESS );
                case IQI_Type::SHEAR_STRESS:
                    return std::make_shared< IQI_Stress >( Stress_Type::SHEAR_STRESS );
                case IQI_Type::VON_MISES_STRESS:
                    return std::make_shared< IQI_Stress >( Stress_Type::VON_MISES_STRESS );
                case IQI_Type::PRINCIPAL_STRESS:
                    return std::make_shared< IQI_Stress >( Stress_Type::PRINCIPAL_STRESS );
                case IQI_Type::STRESS_VECTOR:
                    return std::make_shared< IQI_Stress >( Stress_Type::STRESS_VECTOR );

                case IQI_Type::NORMAL_STRESS_CAUCHY:
                    return std::make_shared< IQI_Stress >( Stress_Type::NORMAL_STRESS, CM_Function_Type::CAUCHY );
                case IQI_Type::SHEAR_STRESS_CAUCHY:
                    return std::make_shared< IQI_Stress >( Stress_Type::SHEAR_STRESS, CM_Function_Type::CAUCHY );
                case IQI_Type::VON_MISES_STRESS_CAUCHY:
                    return std::make_shared< IQI_Stress >( Stress_Type::VON_MISES_STRESS, CM_Function_Type::CAUCHY );
                case IQI_Type::PRINCIPAL_STRESS_CAUCHY:
                    return std::make_shared< IQI_Stress >( Stress_Type::PRINCIPAL_STRESS, CM_Function_Type::CAUCHY );
                case IQI_Type::STRESS_VECTOR_CAUCHY:
                    return std::make_shared< IQI_Stress >( Stress_Type::STRESS_VECTOR, CM_Function_Type::CAUCHY );

                case IQI_Type::MAX_NORMAL_STRESS:
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::NORMAL_STRESS );
                case IQI_Type::MAX_SHEAR_STRESS:
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::SHEAR_STRESS );
                case IQI_Type::MAX_VON_MISES_STRESS:
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::VON_MISES_STRESS );
                case IQI_Type::MAX_PRINCIPAL_STRESS:
                    return std::make_shared< IQI_Max_Stress >( Stress_Type::PRINCIPAL_STRESS );

                case IQI_Type::L2_ERROR_ANALYTIC:
                    return std::make_shared< IQI_L2_Error_Analytic >();

                case IQI_Type::H1_ERROR_ANALYTIC:
                    return std::make_shared< IQI_H1_Error_Analytic >();

                case IQI_Type::H1_ERROR:
                    return std::make_shared< IQI_H1_Error >();

                case IQI_Type::J_INTEGRAL:
                    return std::make_shared< IQI_J_Integral >();

                case IQI_Type::DRAG_COEFF:
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( 1 );

                case IQI_Type::LIFT_COEFF:
                    return std::make_shared< IQI_Drag_Lift_Coefficient >( -1 );

                case IQI_Type::TURBULENT_DYNAMIC_VISCOSITY:
                    return std::make_shared< IQI_Turbulent_Dynamic_Viscosity >();
                case IQI_Type::EFFECTIVE_DYNAMIC_VISCOSITY:
                    return std::make_shared< IQI_Effective_Dynamic_Viscosity >();
                case IQI_Type::EFFECTIVE_CONDUCTIVITY:
                    return std::make_shared< IQI_Effective_Conductivity >();
                case IQI_Type::SPALART_ALLMARAS_COEFFICIENT:
                    return std::make_shared< IQI_Spalart_Allmaras_Coefficient >();

                case IQI_Type::LATENT_HEAT_ABSORPTION:
                    return std::make_shared< IQI_Latent_Heat_Absorption >();

                case IQI_Type::TOTAL_PRESSURE:
                    return std::make_shared< IQI_Total_Pressure >();

                case IQI_Type::MASS_FLOW:
                    return std::make_shared< IQI_Mass_Flow >();

                case IQI_Type::THERMAL_ENERGY_CONVECTIVE_FLUX:
                    return std::make_shared< IQI_Thermal_Energy_Convective_Flux >();
                case IQI_Type::THERMAL_ENERGY_DIFFUSIVE_FLUX:
                    return std::make_shared< IQI_Thermal_Energy_Diffusive_Flux >();
                case IQI_Type::HOMOGENIZED_CONSTITUTIVE:
                    return std::make_shared< IQI_Homogenized_Constitutive >();
                case IQI_Type::POWER_DISSIPATION:
                    return std::make_shared< IQI_Power_Dissipation >();
                case IQI_Type::POWER_DISSIPATION_BULK:
                    return std::make_shared< IQI_Power_Dissipation_Bulk >();

                case IQI_Type::HEAT_METHOD_PENALTY:
                    return std::make_shared< IQI_Heat_Method_Penalty >();

                case IQI_Type::ADVECTION_STRONG_RESIDUAL:
                    return std::make_shared< IQI_Advection_Strong_Residual >();
                case IQI_Type::STRONG_RESIDUAL_SA:
                    return std::make_shared< IQI_Strong_Residual_SA >();
                case IQI_Type::STRONG_RESIDUAL_INCOMPRESSIBLE_NS:
                    return std::make_shared< IQI_Strong_Residual_Incompressible_NS >();

                case IQI_Type::SP_CROSSWIND_INCOMPRESSIBLE_NS:
                    return std::make_shared< IQI_SP_Crosswind_Incompressible_NS >();
                case IQI_Type::RES_CROSSWIND_INCOMPRESSIBLE_NS:
                    return std::make_shared< IQI_Res_Crosswind_Incompressible_NS >();
                case IQI_Type::RES_SUPG_INCOMPRESSIBLE_NS:
                    return std::make_shared< IQI_Res_SUPG_Incompressible_NS >();
                case IQI_Type::SP_CROSSWIND_SA:
                    return std::make_shared< IQI_SP_Crosswind_SA >();

                case IQI_Type::ZIENKIEWICZ_ZHU_VON_MISES_STRESS:
                    return std::make_shared< IQI_Zienkiewicz_Zhu >( Stress_Type::VON_MISES_STRESS );
                case IQI_Type::NORMAL_VECTOR:
                    return std::make_shared< IQI_Normal_Vector >();
                case IQI_Type::JUMP_DOF:
                    return std::make_shared< IQI_Jump_Dof >();
                case IQI_Type::JUMP_TRACTION:
                    return std::make_shared< IQI_Jump_Traction >();
                case IQI_Type::TRACTION:
                    return std::make_shared< IQI_Traction >();
                case IQI_Type::CONTACT_PRESSURE_REFERENCE:
                    return std::make_shared< IQI_Contact_Pressure >(true);
                case IQI_Type::CONTACT_PRESSURE_CURRENT:
                    return std::make_shared< IQI_Contact_Pressure >(false);
                case IQI_Type::GAP:
                    return std::make_shared< IQI_Gap >();

                case IQI_Type::LINEAR_ELASTICITY_DAMAGE:
                    return std::make_shared< IQI_Linear_Elasticity_Damage >();

                default:
                    MORIS_ERROR( false, " IQI_Factory::create_IQI - No IQI type specified. " );
                    return nullptr;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
