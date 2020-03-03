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
            QUAD_1x1,
            QUAD_2x2,
            QUAD_3x3,
            QUAD_4x4,
            QUAD_5x5,
            TRI_1,
            TRI_3,
            TRI_6,
            TRI_7,
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
            SPATIALDIFF_DIRICHLET, // spatial diffusion Dirichlet (Nitsche)
            SPATIALDIFF_NEUMANN,   // spatial diffusion Neumann
            SPATIALDIFF_INTERFACE, // spatial diffusion Nitsche interface condition
            SPATIALDIFF_GHOST,     // spatial diffusion ghost
            SPATIALDIFF_VW_GHOST,  // spatial diffusion virtual work ghost
            STRUC_LINEAR_BULK,     // linear elasticity bulk
            STRUC_LINEAR_DIRICHLET,// linear elasticity Dirichlet (Nitsche)
            STRUC_LINEAR_NEUMANN,  // linear elasticity Neumann
            STRUC_LINEAR_INTERFACE,// linear elasticity Nitsche interface condition
            STRUC_LINEAR_GHOST,    // linear elasticity Ghost field based
            STRUC_LINEAR_VW_GHOST, // linear elasticity Ghost flux based
            STRUC_LINEAR_PRESSURE_BULK, //linear elasticity bulk mixed formulation
            STRUC_LINEAR_PRESSURE_DIRICHLET, // linear elasticity Dirichlet mixed formulation (Nitsche)
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
            PROPERTY,
            L2_ERROR_ANALYTIC,
            H1_ERROR_ANALYTIC,
            H1_SEMI_ERROR,
            J_INTEGRAL,
            END_IQI_TYPE
        };

//------------------------------------------------------------------------------
        enum class Constitutive_Type
        {
            UNDEFINED,
            DIFF_LIN_ISO,
            STRUC_LIN_ISO,
            STRUC_LIN_ISO_PRESSURE,
            FLUID_INCOMPRESSIBLE,
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
            GHOST_DISPL,
            GHOST_VW,
            NITSCHE_INTERFACE,
            MASTER_WEIGHT_INTERFACE,
            SLAVE_WEIGHT_INTERFACE,
            RECIPROCAL_TOTAL_VOLUME,
            END_STABILIZATION_TYPE
        };

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_CL_FEM_ENUMS_HPP_ */
