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

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        enum class Element_Type
        {
            UNDEFINED,
            BULK,
            SIDESET,
            DOUBLE_SIDESET,
            NONCONFORMAL_SIDESET,
            TIME_SIDESET,
            TIME_BOUNDARY,
            END_ELEMENT_TYPE
        };

        //------------------------------------------------------------------------------

        enum class IWG_Type
        {
            UNDEFINED,
            L2,           // L2 projection
            HJ,           // Hamilton-Jacobi
            HJTEST,       // Space time IWG for test only
            HELMHOLTZ,    // Helmholtz
            HELMHOLTZ_INTERFACE_SYMMETRIC_NITSCHE,
            HELMHOLTZ_INTERFACE_UNSYMMETRIC_NITSCHE,
            LSNORMAL,         // LS normal
            OLSSON,           // Olsson et al. (2007) reinitialization
            NONLOCAL_BULK,    // Nonlocal
            HISTORY_BULK,     // History
            NONLOCAL_INTERFACE_SYMMETRIC_NITSCHE,
            NONLOCAL_INTERFACE_UNSYMMETRIC_NITSCHE,

            SPATIALDIFF_BULK,                             // spatial diffusion bulk
            SPATIALDIFF_PC_BULK,                          // spatial diffusion bulk with phase change
            ADVECTION_BULK,                               // spatial advection bulk
            SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE,      // spatial diffusion Dirichlet (Nitsche)
            SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE,    // spatial diffusion Dirichlet (Nitsche)
            SPATIALDIFF_ROBIN_SYMMETRIC_NITSCHE,          // spatial diffusion Robin (Nitsche)
            SPATIALDIFF_ROBIN_UNSYMMETRIC_NITSCHE,        // spatial diffusion Robin (Nitsche)
            SPATIALDIFF_NEUMANN,                          // spatial diffusion Neumann
            SPATIALDIFF_CONVECTION,                       // spatial diffusion Robin (Convection)
            SPATIALDIFF_RADIATION,                        // spatial diffusion Radiation BC
            SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE,      // spatial diffusion Nitsche interface condition
            SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE,    // spatial diffusion Nitsche interface condition
            SPATIALDIFF_GGLS_PC,                          // spatial diffusion GGLS stabilization term for phase change
            SPATIALDIFF_VW_GHOST,                         // spatial diffusion virtual work ghost

            STRUC_LINEAR_BULK,                             // linear elasticity bulk
            STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE,      // linear elasticity Dirichlet (Nitsche)
            STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE,    // linear elasticity Dirichlet (Nitsche)
            STRUC_LINEAR_NEUMANN,                          // linear elasticity Neumann
            STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE,      // linear elasticity Nitsche interface condition
            STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE,    // linear elasticity Nitsche interface condition

            Struc_Linear_Interface_SLM_Constraint,    // linear elasticity stabilzed lagrange multiplier constraint
            Struc_Linear_Interface_SLM_L2,            // linear elasticity stabilzed lagrange multiplier L2 projection
            Struc_Linear_Interface_SLM_Mixed,         // linear elasticity stabilzed lagrange multiplier mixed term
            Struc_Linear_Interface_SLM_LMJump,        // linear elasticity stabilzed lagrange multiplier LM jump term

            STRUC_LINEAR_VW_GHOST,                                  // linear elasticity Ghost flux based
            STRUC_LINEAR_PRESSURE_BULK,                             // linear elasticity bulk mixed formulation
            STRUC_LINEAR_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE,      // linear elasticity Dirichlet mixed formulation (Nitsche)
            STRUC_LINEAR_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE,    // linear elasticity Dirichlet mixed formulation (Nitsche)
            STRUC_LINEAR_FLUID_INTERFACE,                           // one-side fluid-structure coupling
            STRUC_VON_MISES_STRESS,

            STRUC_NON_LINEAR_BULK_SE,                                    // nonlinear elasticity bulk, based on S and E work conjugates
            STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_SE,             // nonlinear elasticity sym Dirichlet (Nitsche), based on S and E work conjugates
            STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_SE,           // nonlinear elasticity nonsym Dirichlet (Nitsche), based on S and E work conjugates
            STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_SE,             // nonlinear elasticity sym Interface (Nitsche), based on S and E work conjugates
            STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_SE,           // nonlinear elasticity nonsym Interface (Nitsche), based on S and E work conjugates
            STRUC_NON_LINEAR_BULK_PF,                                    // nonlinear elasticity bulk, based on P and F work conjugates
            STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF,             // nonlinear elasticity sym Dirichlet (Nitsche), based on P and F work conjugates
            STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_PF,           // nonlinear elasticity nonsym Dirichlet (Nitsche), based on P and F work conjugates
            STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_PF,             // nonlinear elasticity sym Interface (Nitsche), based on P and F work conjugates
            STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_PF,           // nonlinear elasticity nonsym Interface (Nitsche), based on P and F work conjugates
            STRUC_NON_LINEAR_BULK_CAUCHYEPS,                             // nonlinear elasticity bulk, based on Cauchy and Epsilon work conjugates
            STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_CAUCHYEPS,      // nonlinear elasticity sym Dirichlet (Nitsche), based on Cauchy and Epsilon work conjugates
            STRUC_NON_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE_CAUCHYEPS,    // nonlinear elasticity nonsym Dirichlet (Nitsche), based on Cauchy and Epsilon work conjugates
            STRUC_NON_LINEAR_INTERFACE_SYMMETRIC_NITSCHE_CAUCHYEPS,      // nonlinear elasticity sym Interface (Nitsche), based on Cauchy and Epsilon work conjugates
            STRUC_NON_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE_CAUCHYEPS,    // nonlinear elasticity nonsym Interface (Nitsche), based on Cauchy and Epsilon work conjugates

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
            STRUC_LINEAR_CONTACT_GAP_SYMMETRIC_NITSCHE_UNBIASED,
            STRUC_LINEAR_CONTACT_GAP_UNSYMMETRIC_NITSCHE_UNBIASED,
            STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED,
            STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED,
            STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED,
            STRUC_LINEAR_CONTACT_PENALTY,

            STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_SYMMETRIC,
            STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_UNSYMMETRIC,
            STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_NEUTRAL,
            STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_SYMMETRIC,
            STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_UNSYMMETRIC,
            STRUC_NONLINEAR_CONTACT_SEITZ_UNBIASED_NEUTRAL,


            GHOST_NORMAL_FIELD,

            DEBUG,
            END_IWG_TYPE
        };

        //------------------------------------------------------------------------------

        enum class IQI_Type
        {
            UNDEFINED,
            VOLUME,    // volume
            STRAIN_ENERGY,
            VOLUME_FRACTION,
            RAY_LENGTH,
            DOF,
            EIGEN_VECTOR,
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
            NORMAL_VECTOR,

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

            END_IQI_TYPE
        };
        //------------------------------------------------------------------------------

        enum class Constitutive_Type
        {
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
            STRUC_NON_LIN_ISO_NEO_HOOKEAN,
            FLUID_INCOMPRESSIBLE,
            FLUID_TURBULENCE,
            FLUID_COMPRESSIBLE_IDEAL,
            FLUID_COMPRESSIBLE_VDW,
            FLUID_COMPRESSIBLE_NEWTONIAN,
            SPALART_ALLMARAS_TURBULENCE,
            END_CONSTITUTIVE_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Material_Type
        {
            UNDEFINED,
            PERFECT_GAS,
            VAN_DER_WAALS_FLUID,
            END_MATERIAL_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Variable_Set
        {
            UNDEFINED,
            CONSERVATIVE,
            DENSITY_PRIMITIVE,
            PRESSURE_PRIMITIVE,
            ENTROPY,
            END_VARIABLE_SET
        };

        //------------------------------------------------------------------------------

        enum class Model_Type
        {
            UNDEFINED,
            PLANE_STRESS,
            PLANE_STRAIN,
            AXISYMMETRIC,
            FULL,
            HYDROSTATIC,    // not implemented yet
            DEVIATORIC,
            END_MODEL_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Stabilization_Type
        {
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
            END_STABILIZATION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Measure_Type
        {
            UNDEFINED,
            CELL_MEASURE,
            CELL_SIDE_MEASURE,
            CELL_LENGTH_MEASURE,
            END_MEASURE_TYPE
        };

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

        enum class FDScheme_Type
        {
            UNDEFINED,
            POINT_1_FORWARD,
            POINT_1_BACKWARD,
            POINT_3_CENTRAL,
            POINT_5,
            END_FD_SCHEME
        };

        //------------------------------------------------------------------------------

        enum class Perturbation_Type
        {
            UNDEFINED,
            RELATIVE,
            ABSOLUTE,
            END_PERTURBATION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Stress_Type
        {
            UNDEFINED,
            NORMAL_STRESS,
            SHEAR_STRESS,
            VON_MISES_STRESS,
            PRINCIPAL_STRESS,
            MAX_SHEAR_STRESS,
            STRESS_VECTOR,
            END_STRESS_TYPE
        };

        //------------------------------------------------------------------------------

        enum class CM_Function_Type
        {
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
            LAGRANGIAN,
            EULERIAN,
            END_CM_FUNCTION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Time_Continuity_Flag
        {
            DEFAULT,
            TIME_CONTINUITY_ONLY,
            NO_TIME_CONTINUITY,
            END_TIME_CONTINUITY_FLAG,
        };

        //------------------------------------------------------------------------------

        enum class CM_Request_Type
        {
            UNDEFINED,
            STRAIN,
            FLUX,
            TRACTION,
            TEST_TRACTION,
            DAMAGE,
            SMOOTH_DAMAGE,
            EQSTRAIN,
            HISTORY,
            END_CM_REQUEST_TYPE
        };
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ENUMS_HPP_ */
