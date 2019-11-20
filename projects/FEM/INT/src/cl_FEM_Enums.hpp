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
            BEZIER    // Bezier type elements
        };

//------------------------------------------------------------------------------

        enum class Integration_Type
        {
            UNDEFINED,
            CONSTANT,
            GAUSS    // Gauss ( Quad and Hex ), Dunavant ( Tri ), Hammer ( Tet )
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
            ST_QUAD_3x3x1 // < -- space time combined
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
            HJTEST,     // FIXME for test only, to be deleted
            HELMHOLTZ,  // Helmholtz
            LSNORMAL,   // LS normal
            OLSSON,     // Olsson et al. (2007) reinitialization
            SPATIALDIFF_BULK,      // spatial diffusion bulk
            SPATIALDIFF_DIRICHLET, // spatial diffusion Dirichlet (Nitsche)
            SPATIALDIFF_NEUMANN,   // spatial diffusion Neumann
            SPATIALDIFF_GHOST,     // spatial diffusion ghost
            SPATIALDIFF_VW_GHOST,  // spatial diffusion virtual work ghost
            SPATIALDIFF_INTERFACE, // spatial diffusion Nitsche interface condition
            STRUC_LINEAR_BULK,     // linear elasticity bulk
            STRUC_LINEAR_DIRICHLET,// linear elasticity Dirichlet (Nitsche)
            STRUC_LINEAR_NEUMANN,  // linear elasticity Neumann
            STRUC_LINEAR_INTERFACE,// linear elasticity Nitsche interface condition
            END_IWG_TYPE
        };

//------------------------------------------------------------------------------

        enum class Property_Type
        {
            UNDEFINED,
            TEMP_DIRICHLET,
            TEMP_NEUMANN,
            CONDUCTIVITY,
            TEMP_LOAD,
            STRUC_DIRICHLET,
            STRUC_NEUMANN,
            YOUNGS_MODULUS,
            POISSONS_RATIO,
            END_PROPERTY_TYPE
        };

//------------------------------------------------------------------------------

        enum class Constitutive_Type
        {
            UNDEFINED,
            DIFF_LIN_ISO,
            STRUC_LIN_ISO,
            END_CONSTITUTIVE_TYPE
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
    } /* namespace fem */
} /* namespace moris */


#endif /* SRC_FEM_CL_FEM_ENUMS_HPP_ */
