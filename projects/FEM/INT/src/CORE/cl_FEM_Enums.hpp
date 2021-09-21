/*
 * cl_FEM_Enums.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
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
                TIME_SIDESET,
                TIME_BOUNDARY,
                END_ELEMENT_TYPE
        };

        //------------------------------------------------------------------------------

        enum class IWG_Type
        {
                UNDEFINED,
                L2,         // L2 projection
                HJ,         // Hamilton-Jacobi
                HJTEST,     // Space time IWG for test only
                HELMHOLTZ,  // Helmholtz
                LSNORMAL,   // LS normal
                OLSSON,     // Olsson et al. (2007) reinitialization

                SPATIALDIFF_BULK,      // spatial diffusion bulk
                SPATIALDIFF_PC_BULK,   // spatial diffusion bulk with phase change
                ADVECTION_BULK,
                SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE, // spatial diffusion Dirichlet (Nitsche)
                SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE,
                SPATIALDIFF_NEUMANN,   // spatial diffusion Neumann
                SPATIALDIFF_ROBIN,   // spatial diffusion Robin (Convection)
                SPATIALDIFF_RADIATION,   // spatial diffusion Radiation BC
                SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE, // spatial diffusion Nitsche interface condition
                SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE,
                SPATIALDIFF_GGLS_PC,   // spatial diffusion GGLS stabilization term for phase change
                SPATIALDIFF_VW_GHOST,  // spatial diffusion virtual work ghost

                STRUC_LINEAR_BULK,     // linear elasticity bulk
                STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE,// linear elasticity Dirichlet (Nitsche)
                STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE,// linear elasticity Dirichlet (Nitsche)
                STRUC_LINEAR_NEUMANN,  // linear elasticity Neumann
                STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE,// linear elasticity Nitsche interface condition
                STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE,
                STRUC_LINEAR_VW_GHOST, // linear elasticity Ghost flux based
                STRUC_LINEAR_PRESSURE_BULK, //linear elasticity bulk mixed formulation
                STRUC_LINEAR_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE, // linear elasticity Dirichlet mixed formulation (Nitsche)
                STRUC_LINEAR_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                STRUC_VON_MISES_STRESS,

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
                STRUC_LINEAR_CONTACT_PENALTY,
                GHOST_NORMAL_FIELD,
                END_IWG_TYPE
        };

        //------------------------------------------------------------------------------

        enum class IQI_Type
        {
                UNDEFINED,
                VOLUME,         // volume
                STRAIN_ENERGY,
                VOLUME_FRACTION,
                DOF,
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
                TURBULENT_KINEMATIC_VISCOSITY,
                TOTAL_PRESSURE,
                MASS_FLOW,
                THERMAL_ENERGY_CONVECTIVE_FLUX,
                THERMAL_ENERGY_DIFFUSIVE_FLUX,

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

                HOMOGENIZED_CONSTITUTIVE,
                HEAT_METHOD_PENALTY,

                END_IQI_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Constitutive_Type
        {
                UNDEFINED,
                DIFF_LIN_ISO,
                DIFF_LIN_ISO_PC,
                STRUC_LIN_ISO,
                STRUC_LIN_ISO_PRESSURE,
                STRUC_NONLIN_ISO,
                FLUID_INCOMPRESSIBLE,
                FLUID_TURBULENCE,
                FLUID_COMPRESSIBLE_IDEAL,
                FLUID_COMPRESSIBLE_VDW,
                FLUID_COMPRESSIBLE_NEWTONIAN,
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
                HYDROSTATIC, // not implemented yet
                DEVIATORIC,
                END_MODEL_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Stabilization_Type
        {
            UNDEFINED,
            DIRICHLET_NITSCHE,
            GGLS_DIFFUSION,
            GHOST_DISPL,
            GHOST_NORMAL_FIELD,
            GHOST_VW,
            NITSCHE_INTERFACE,
            MASTER_WEIGHT_INTERFACE,
            SLAVE_WEIGHT_INTERFACE,
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
            TURBULENCE_DIRICHLET_NITSCHE,
            SPALART_ALLMARAS_NITSCHE_INTERFACE,
            PENALTY_CONTACT,
            STAB_PENALTY_CONTACT,
            MEASURE,
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

        inline
        map< std::string, enum fem::Measure_Type > get_measure_type_map()
        {
            map< std::string, enum fem::Measure_Type > tFemMeasureTypeMap;

            tFemMeasureTypeMap["UNDEFINED"]           = fem::Measure_Type::UNDEFINED;
            tFemMeasureTypeMap["CELL_MEASURE"]        = fem::Measure_Type::CELL_MEASURE;
            tFemMeasureTypeMap["CELL_SIDE_MEASURE"]   = fem::Measure_Type::CELL_SIDE_MEASURE;
            tFemMeasureTypeMap["CELL_LENGTH_MEASURE"] = fem::Measure_Type::CELL_LENGTH_MEASURE;
            tFemMeasureTypeMap["END_MEASURE_TYPE"]    = fem::Measure_Type::END_MEASURE_TYPE;

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
                END_CM_FUNCTION_TYPE
        };

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ENUMS_HPP_ */
