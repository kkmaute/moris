/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_Enums.hpp
 *
 */

#ifndef SRC_FEM_CL_FEM_ENUMS_HPP_
#define SRC_FEM_CL_FEM_ENUMS_HPP_

#include "assert.hpp"
#include "cl_Map.hpp"
#include "fn_enum_macros.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        ENUM_MACRO( Element_Type,
                UNDEFINED,
                BULK,
                SIDESET,
                DOUBLE_SIDESET,
                TIME_SIDESET,
                TIME_BOUNDARY,
                END_ELEMENT_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( IWG_Type,
                UNDEFINED,
                L2,
                HJ,
                HJTEST,
                HELMHOLTZ,
                HELMHOLTZ_INTERFACE_SYMMETRIC_NITSCHE,
                HELMHOLTZ_INTERFACE_UNSYMMETRIC_NITSCHE,
                LSNORMAL,
                OLSSON,
                NONLOCAL_INTERFACE_SYMMETRIC_NITSCHE,
                NONLOCAL_INTERFACE_UNSYMMETRIC_NITSCHE,
                L2_EQSTRAIN_BULK,
                L2_HISTORY_BULK,
                L2_DAMAGE_BULK,
                SPATIALDIFF_BULK,
                SPATIALDIFF_PC_BULK,
                ADVECTION_BULK,
                SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE,
                SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE,
                SPATIALDIFF_ROBIN_SYMMETRIC_NITSCHE,
                SPATIALDIFF_ROBIN_UNSYMMETRIC_NITSCHE,
                SPATIALDIFF_NEUMANN,
                SPATIALDIFF_CONVECTION,
                SPATIALDIFF_RADIATION,
                SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE,
                SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE,
                SPATIALDIFF_GGLS_PC,
                SPATIALDIFF_VW_GHOST,
                STRUC_LINEAR_BULK,
                STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE,
                STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_NEUMANN,
                STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE,
                STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE,
                Struc_Linear_Interface_SLM_Constraint,
                Struc_Linear_Interface_SLM_L2,
                Struc_Linear_Interface_SLM_Mixed,
                Struc_Linear_Interface_SLM_LMJump,
                STRUC_LINEAR_VW_GHOST,
                STRUC_LINEAR_PRESSURE_BULK,
                STRUC_LINEAR_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE,
                STRUC_LINEAR_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_FLUID_INTERFACE,
                STRUC_VON_MISES_STRESS,
                STRUC_NON_LINEAR_BULK_SE,
                STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS,
                STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_SE,
                STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_SE,
                STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_SE,
                STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_SE,
                STRUC_NON_LINEAR_BULK_PF,
                STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF,
                STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF,
                STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_PF,
                STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_PF,
                STRUC_NON_LINEAR_BULK_CAUCHYEPS,
                STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_CAUCHYEPS,
                STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_CAUCHYEPS,
                STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_CAUCHYEPS,
                STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_CAUCHYEPS,
                INCOMPRESSIBLE_NS_VELOCITY_BULK,
                INCOMPRESSIBLE_NS_PRESSURE_BULK,
                INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST,
                INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_SLIPBOUNDARY_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_SLIPBOUNDARY_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_SLIPBOUNDARY_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_SLIPBOUNDARY_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_IMPOSED_PRESSURE,
                INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_UNSYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_BULK,
                COMPRESSIBLE_NS_BOUNDARY,
                COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_DENSITY_BULK,
                COMPRESSIBLE_NS_VELOCITY_BULK,
                COMPRESSIBLE_NS_TEMPERATURE_BULK,
                COMPRESSIBLE_NS_ADVECTIVE_MOMENTUM_FLUX,
                COMPRESSIBLE_NS_ADVECTIVE_ENERGY_FLUX,
                COMPRESSIBLE_NS_MASS_FLUX_NEUMANN,
                COMPRESSIBLE_NS_TRACTION_NEUMANN,
                COMPRESSIBLE_NS_HEAT_FLUX_NEUMANN,
                COMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_TEMPERATURE_DIRICHLET_SYMMETRIC_NITSCHE,
                COMPRESSIBLE_NS_TEMPERATURE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                FS_STRUC_INTERFACE,
                TIME_CONTINUITY_DOF,
                SPALART_ALLMARAS_TURBULENCE_BULK,
                SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE,
                SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                SPALART_ALLMARAS_TURBULENCE_INTERFACE_SYMMETRIC_NITSCHE,
                SPALART_ALLMARAS_TURBULENCE_INTERFACE_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_CONTACT_SYMMETRIC_NITSCHE,
                STRUC_LINEAR_CONTACT_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_CONTACT_GAP_SYMMETRIC_NITSCHE,
                STRUC_LINEAR_CONTACT_GAP_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_CONTACT_PENALTY,
                GHOST_NORMAL_FIELD,
                USER_DEFINED,
                END_IWG_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( IQI_Type,
                UNDEFINED,
                VOLUME,
                STRAIN_ENERGY,
                VOLUME_FRACTION,
                DOF,
                EIGEN_VECTOR,
                EIGEN_VALUE,
                ALM_DOF,
                MAX_DOF,
                PROPERTY,
                STABILIZATION,
                L2_ERROR_ANALYTIC,
                H1_ERROR_ANALYTIC,
                H1_ERROR,
                J_INTEGRAL,
                LIFT_COEFF,
                DRAG_COEFF,
                LATENT_HEAT_ABSORPTION,
                TURBULENT_DYNAMIC_VISCOSITY,
                EFFECTIVE_DYNAMIC_VISCOSITY,
                EFFECTIVE_CONDUCTIVITY,
                SPALART_ALLMARAS_COEFFICIENT,
                POWER_DISSIPATION,
                POWER_DISSIPATION_BULK,
                TOTAL_PRESSURE,
                MASS_FLOW,
                THERMAL_ENERGY_CONVECTIVE_FLUX,
                THERMAL_ENERGY_DIFFUSIVE_FLUX,
                JUMP_DOF,
                JUMP_TRACTION,
                TRACTION,
                ADVECTION_STRONG_RESIDUAL,
                STRONG_RESIDUAL_SA,
                STRONG_RESIDUAL_INCOMPRESSIBLE_NS,
                SP_CROSSWIND_INCOMPRESSIBLE_NS,
                SP_CROSSWIND_SA,
                RES_CROSSWIND_INCOMPRESSIBLE_NS,
                RES_SUPG_INCOMPRESSIBLE_NS,
                MAX_STRESS,
                MAX_NORMAL_STRESS,
                MAX_SHEAR_STRESS,
                MAX_VON_MISES_STRESS,
                MAX_PRINCIPAL_STRESS,
                STRESS,
                NORMAL_STRESS,
                SHEAR_STRESS,
                VON_MISES_STRESS,
                PRINCIPAL_STRESS,
                STRESS_VECTOR,
                NORMAL_STRESS_CAUCHY,
                SHEAR_STRESS_CAUCHY,
                VON_MISES_STRESS_CAUCHY,
                PRINCIPAL_STRESS_CAUCHY,
                STRESS_VECTOR_CAUCHY,
                HOMOGENIZED_CONSTITUTIVE,
                HEAT_METHOD_PENALTY,
                MIN_FEATURE_SIZE,
                ZIENKIEWICZ_ZHU_VON_MISES_STRESS,
                LINEAR_ELASTICITY_DAMAGE,
                END_IQI_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Constitutive_Type,
                UNDEFINED,
                DIFF_LIN_ISO,
                DIFF_LIN_ISO_PC,
                DIFF_LIN_ISO_TURBULENCE,
                STRUC_LIN_ISO,
                STRUC_LIN_MT,
                STRUC_LIN_ISO_PRESSURE,
                STRUC_LIN_ISO_DAMAGE,
                STRUC_NON_LIN_ISO,
                STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF,
                STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_BONET,
                STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS,
                FLUID_INCOMPRESSIBLE,
                FLUID_TURBULENCE,
                FLUID_COMPRESSIBLE_IDEAL,
                FLUID_COMPRESSIBLE_VDW,
                FLUID_COMPRESSIBLE_NEWTONIAN,
                SPALART_ALLMARAS_TURBULENCE,
                END_CONSTITUTIVE_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Material_Type,
                UNDEFINED,
                PERFECT_GAS,
                VAN_DER_WAALS_FLUID,
                END_MATERIAL_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Variable_Set,
                UNDEFINED,
                CONSERVATIVE,
                DENSITY_PRIMITIVE,
                PRESSURE_PRIMITIVE,
                ENTROPY,
                END_VARIABLE_SET )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Model_Type,
                UNDEFINED,
                PLANE_STRESS,
                PLANE_STRAIN,
                AXISYMMETRIC,
                FULL,
                HYDROSTATIC,    // not implemented yet
                DEVIATORIC,
                END_MODEL_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Stabilization_Type,
                UNDEFINED,
                DIRICHLET_NITSCHE,
                ROBIN_NITSCHE,
                GGLS_DIFFUSION,
                GHOST_DISPL,
                GHOST_NORMAL_FIELD,
                GHOST_VW,
                NITSCHE_INTERFACE,
                LEADER_WEIGHT_INTERFACE,
                FOLLOWER_WEIGHT_INTERFACE,
                RECIPROCAL_TOTAL_VOLUME,
                INCOMPRESSIBLE_FLOW,
                VISCOUS_GHOST,
                CONVECTIVE_GHOST,
                PRESSURE_GHOST,
                TIME_VELOCITY_GHOST,
                VELOCITY_DIRICHLET_NITSCHE,
                VELOCITY_SLIPBOUNDARY_NITSCHE,
                COMPRESSIBLE_VELOCITY_DIRICHLET_NITSCHE,
                COMPRESSIBLE_DIRICHLET_NITSCHE,
                SUPG_ADVECTION,
                SUPG_SPALART_ALLMARAS_TURBULENCE,
                YZBETA_ADVECTION,
                CROSSWIND,
                TURBULENCE_DIRICHLET_NITSCHE,
                SPALART_ALLMARAS_NITSCHE_INTERFACE,
                PENALTY_CONTACT,
                STAB_PENALTY_CONTACT,
                MEASURE,
                LAGRANGE_MULTIPLIER_L2,
                END_STABILIZATION_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Measure_Type,
                UNDEFINED,
                CELL_MEASURE,
                CELL_SIDE_MEASURE,
                CELL_LENGTH_MEASURE,
                END_MEASURE_TYPE )

        inline map< std::string, enum fem::Measure_Type >
        get_measure_type_map()
        {
            map< std::string, enum fem::Measure_Type > tFemMeasureTypeMap;

            tFemMeasureTypeMap[ "UNDEFINED" ]           = fem::Measure_Type::UNDEFINED;
            tFemMeasureTypeMap[ "CELL_MEASURE" ]        = fem::Measure_Type::CELL_MEASURE;
            tFemMeasureTypeMap[ "CELL_SIDE_MEASURE" ]   = fem::Measure_Type::CELL_SIDE_MEASURE;
            tFemMeasureTypeMap[ "CELL_LENGTH_MEASURE" ] = fem::Measure_Type::CELL_LENGTH_MEASURE;
            tFemMeasureTypeMap[ "END_MEASURE_TYPE" ]    = fem::Measure_Type::END_MEASURE_TYPE;

            return tFemMeasureTypeMap;
        }

        //------------------------------------------------------------------------------

        ENUM_MACRO( FDScheme_Type,
                UNDEFINED,
                POINT_1_FORWARD,
                POINT_1_BACKWARD,
                POINT_3_CENTRAL,
                POINT_5,
                END_FD_SCHEME )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Perturbation_Type,
                UNDEFINED,
                RELATIVE,
                ABSOLUTE,
                END_PERTURBATION_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Stress_Type,
                UNDEFINED,
                NORMAL_STRESS,
                SHEAR_STRESS,
                VON_MISES_STRESS,
                PRINCIPAL_STRESS,
                MAX_SHEAR_STRESS,
                STRESS_VECTOR,
                END_STRESS_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( CM_Function_Type,
                DEFAULT,
                THERMAL,
                FLUID,
                MECHANICAL,
                ENERGY,
                WORK,
                HEAT,
                PRESSURE,
                PK1,
                PK2,
                CAUCHY,
                DEFORMATION_GRADIENT,
                RIGHT_CAUCHY_GREEN,
                INV_RIGHT_CAUCHY_GREEN,
                LEFT_CAUCHY_GREEN,
                INV_LEFT_CAUCHY_GREEN,
                LAGRANGIAN,
                EULERIAN,
                INFINITESIMAL,
                END_CM_FUNCTION_TYPE )

        //------------------------------------------------------------------------------

        ENUM_MACRO( Time_Continuity_Flag,
                DEFAULT,
                TIME_CONTINUITY_ONLY,
                NO_TIME_CONTINUITY,
                GEOMETRIC_STIFFNESS_ONLY,
                END_TIME_CONTINUITY_FLAG )

        //------------------------------------------------------------------------------

        ENUM_MACRO( CM_Request_Type,
                UNDEFINED,
                STRAIN,
                TEST_STRAIN,
                FLUX,
                TRACTION,
                TEST_TRACTION,
                DAMAGE,
                SMOOTH_DAMAGE,
                EQSTRAIN,
                HISTORY,
                END_CM_REQUEST_TYPE )
    }    // namespace fem
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ENUMS_HPP_ */
