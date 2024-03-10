/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Factory.cpp
 *
 */

#include "assert.hpp"
//FEM/INT/src
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Phase_Change.hpp"
#include "cl_FEM_CM_Diffusion_Linear_Isotropic_Turbulence.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"
#include "cl_FEM_CM_Struc_Linear_Isotropic_Damage.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Bonet.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Wriggers.hpp"
#include "cl_FEM_CM_Fluid_Incompressible.hpp"
#include "cl_FEM_CM_Fluid_Turbulence.hpp"
#include "cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp"
#include "cl_FEM_CM_Fluid_Compressible_Ideal.hpp"
#include "cl_FEM_CM_Compressible_Newtonian_Fluid.hpp"
#include "cl_FEM_CM_Fluid_Compressible_Van_der_Waals.hpp"
#include "cl_FEM_CM_Struc_Linear_MoriTanaka.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        std::shared_ptr< Constitutive_Model > CM_Factory::create_CM( fem::Constitutive_Type aConstitutiveType )
        {

            switch( aConstitutiveType )
            {
                case Constitutive_Type::DIFF_LIN_ISO :
                    return std::make_shared< CM_Diffusion_Linear_Isotropic >();

                case Constitutive_Type::DIFF_LIN_ISO_PC :
                    return std::make_shared< CM_Diffusion_Linear_Isotropic_Phase_Change >();

                case Constitutive_Type::DIFF_LIN_ISO_TURBULENCE :
                    return std::make_shared< CM_Diffusion_Linear_Isotropic_Turbulence >();

                case  Constitutive_Type::STRUC_LIN_ISO :
                    return std::make_shared< CM_Struc_Linear_Isotropic >();

                case Constitutive_Type::STRUC_LIN_ISO_DAMAGE:
                    return std::make_shared< CM_Struc_Linear_Isotropic_Damage >();

                case Constitutive_Type::STRUC_LIN_MT:
                    return std::make_shared< CM_Struc_Linear_MoriTanaka >();

                case  Constitutive_Type::STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF :
                    return std::make_shared< CM_Struc_Nonlinear_Isotropic_Saint_Venant_Kirchhoff >();

                case Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_BONET :
                    return std::make_shared< CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Bonet >();

                case Constitutive_Type::STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_WRIGGERS :
                    return std::make_shared< CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Wriggers >();

                case Constitutive_Type::FLUID_INCOMPRESSIBLE :
                    return std::make_shared< CM_Fluid_Incompressible >();

                case Constitutive_Type::FLUID_TURBULENCE :
                    return std::make_shared< CM_Fluid_Turbulence >();

                case Constitutive_Type::SPALART_ALLMARAS_TURBULENCE :
                    return std::make_shared< CM_Spalart_Allmaras_Turbulence >();

                case Constitutive_Type::FLUID_COMPRESSIBLE_IDEAL :
                    return std::make_shared< CM_Fluid_Compressible_Ideal >();

                case Constitutive_Type::FLUID_COMPRESSIBLE_NEWTONIAN :
                    return std::make_shared< CM_Compressible_Newtonian_Fluid >();

                case Constitutive_Type::FLUID_COMPRESSIBLE_VDW :
                    return std::make_shared< CM_Fluid_Compressible_Van_der_Waals >();

                default:
                    MORIS_ERROR( false, " CM_Factory::create_CM - No constitutive type specified. " );
                    return nullptr;
            }
        }
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

