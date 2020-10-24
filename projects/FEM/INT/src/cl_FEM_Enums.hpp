/*
 * cl_FEM_Enums.hpp
 *
 *  Created on: Jul 9, 2018
 *      Author: messe
 */

#ifndef SRC_FEM_CL_FEM_ENUMS_HPP_
#define SRC_FEM_CL_FEM_ENUMS_HPP_

#include "assert.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        enum class Interpolation_Type
        {
                UNDEFINED,
                CONSTANT, // constant interpolation
                LAGRANGE, // the most common finite element types
                BEZIER,   // Bezier type elements
                END_INTERPOLATION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Integration_Type
        {
                UNDEFINED,
                CONSTANT,
                GAUSS,    // Gauss ( Quad and Hex ), Dunavant ( Tri ), Hammer ( Tet )
                END_INTEGRATION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Integration_Order
        {
                UNDEFINED,
                POINT,
                BAR_1,
                BAR_2,
                BAR_3,
                BAR_4,
                BAR_5,
                BAR_6,
                QUAD_1x1,
                QUAD_2x2,
                QUAD_3x3,
                QUAD_4x4,
                QUAD_5x5,
                TRI_1,
                TRI_3,
                TRI_4,
                TRI_6,
                TRI_7,
                TRI_12,
                TRI_13,
                TRI_16,
                TRI_19,
                TRI_25,
                HEX_1x1x1,
                HEX_2x2x2,
                HEX_3x3x3,
                HEX_4x4x4,
                HEX_5x5x5,
                TET_1,
                TET_4,
                TET_5,
                TET_10,
                TET_11,
                TET_15,
                TET_20,
                TET_35,
                TET_56,
                END_INTEGRATION_ORDER
        };

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
                SPATIALDIFF_INTERFACE, // spatial diffusion Nitsche interface condition
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

                INCOMPRESSIBLE_NS_VELOCITY_BULK,
                INCOMPRESSIBLE_NS_PRESSURE_BULK,
                INCOMPRESSIBLE_NS_CONVECTIVE_VELOCITY_GHOST,
                INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_IMPOSED_PRESSURE,
                INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_VELOCITY_INTERFACE_UNSYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_SYMMETRIC_NITSCHE,
                INCOMPRESSIBLE_NS_PRESSURE_INTERFACE_UNSYMMETRIC_NITSCHE,

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
                L2_ERROR_ANALYTIC,
                H1_ERROR_ANALYTIC,
                H1_SEMI_ERROR,
                J_INTEGRAL,
                LIFT_COEFF,
                DRAG_COEFF,
                LATENT_HEAT_ABSORPTION,
                TURBULENT_KINEMATIC_VISCOSITY,
                TOTAL_PRESSURE,
                MASS_FLOW,
                THERMAL_ENERGY,

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

                END_IQI_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Constitutive_Type
        {
                UNDEFINED,
                DIFF_LIN_ISO,
                DIFF_LIN_ISO_PC, // DIFF_LIN_ISO with phase change
                STRUC_LIN_ISO,
                STRUC_LIN_ISO_PRESSURE,
                FLUID_INCOMPRESSIBLE,
                FLUID_TURBULENCE,
                FLUID_COMPRESSIBLE_IDEAL,
                FLUID_COMPRESSIBLE_VDW,
                END_CONSTITUTIVE_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Model_Type
        {
                UNDEFINED,
                PLANE_STRESS,
                PLANE_STRAIN,
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
            SUPG_ADVECTION,
            SUPG_SPALART_ALLMARAS_TURBULENCE,
            TURBULENCE_DIRICHLET_NITSCHE,
            SPALART_ALLMARAS_NITSCHE_INTERFACE,
            PENALTY_CONTACT,
            STAB_PENALTY_CONTACT,
            END_STABILIZATION_TYPE
        };

        //------------------------------------------------------------------------------

        enum class Cluster_Measure
        {
                UNDEFINED,
                MASTER_VOLUME,
                SLAVE_VOLUME,
                INTERFACE_SURFACE,
                ELEMENT_SIZE,
                END_CLUSTER_MEASURE
        };

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
