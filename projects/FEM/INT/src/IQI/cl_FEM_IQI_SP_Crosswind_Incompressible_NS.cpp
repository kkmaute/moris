/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strong_Residual_Incompressible_NS.cpp
 *
 */

#include "cl_FEM_IQI_SP_Crosswind_Incompressible_NS.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_SP_Crosswind_Incompressible_NS::IQI_SP_Crosswind_Incompressible_NS()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::SP_CROSSWIND_INCOMPRESSIBLE_NS;
            init_property( "Gravity", IQI_Property_Type::GRAVITY );
            init_property( "ThermalExpansion", IQI_Property_Type::THERMAL_EXPANSION );
            init_property( "ReferenceTemp", IQI_Property_Type::REF_TEMP );
            init_property( "InvPermeability", IQI_Property_Type::INV_PERMEABILITY );
            init_property( "MassSource", IQI_Property_Type::MASS_SOURCE );
            init_property( "Load", IQI_Property_Type::BODY_LOAD );
            init_constitutive_model( "IncompressibleFluid", IQI_Constitutive_Type::INCOMPRESSIBLE_FLUID );
            init_stabilization_parameter( "Crosswind", IQI_Stabilization_Type::CROSSWIND );
        }

        //------------------------------------------------------------------------------

        void
        IQI_SP_Crosswind_Incompressible_NS::compute_QI( Matrix< DDRMat >& aQI )
        {
            // check if index was set
            if ( mQuantityDofType.size() > 1 )
            {
                MORIS_ERROR( mIQITypeIndex != -1,
                        "IQI_SP_Crosswind_Incompressible_NS::compute_QI - mIQITypeIndex not set." );
            }
            else
            {
                mIQITypeIndex = 0;
            }

            // check if IQI should be computed
            if ( mIQITypeIndex < 0 )
            {
                aQI = MORIS_REAL_MAX;

                return;
            }

            // get the velocity and pressure FIs
            Field_Interpolator* tVelocityFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get space dimension
            uint tSpaceDim = tVelocityFI->get_number_of_fields();

            // check for proper value of IQITypeIndex
            MORIS_ASSERT( (uint)mIQITypeIndex < tSpaceDim,
                    "IQI_SP_Crosswind_Incompressible_NS::compute_QI - incorrect vectorial index." );

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM = get_leader_constitutive_model(IQI_Constitutive_Type::INCOMPRESSIBLE_FLUID);

            // get dynamic viscosity property from CM
            const std::shared_ptr< Property >& tDynViscosityProp =    //
                    tIncFluidCM->get_property( "Viscosity" );

            // get the crosswind stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPCrosswind = get_stabilization_parameter(IQI_Stabilization_Type::CROSSWIND);

            // get zero tolerance for crosswind
            const real tEpsilon = tSPCrosswind->val()( 1 );

            // compute the residual strong form
            Matrix< DDRMat > tRM;
            this->compute_residual_strong_form_momentum( tRM );

            // compute the norm of the gradient of v in direction
            real tNormDir = std::max( norm( tVelocityFI->gradx( 1 )( { 0, tSpaceDim - 1 }, { mIQITypeIndex, mIQITypeIndex } ) ), tEpsilon );

            // compute the abs of the strong form of the residual
            real tRAbsDir = std::max( std::abs( tRM( mIQITypeIndex ) ), tEpsilon );

            // compute full crosswind stabilization parameter value
            aQI = std::max( tSPCrosswind->val()( 0 ) * tRAbsDir / tNormDir - tDynViscosityProp->val()( 0 ), 0.0 );
        }

        //------------------------------------------------------------------------------

        void
        IQI_SP_Crosswind_Incompressible_NS::compute_QI( real aWStar )
        {
            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( get_name() );

            // evaluate strong form
            Matrix< DDRMat > tQI( 1, 1 );
            this->compute_QI( tQI );

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * tQI;
        }

        //------------------------------------------------------------------------------

        void
        IQI_SP_Crosswind_Incompressible_NS::compute_dQIdu( real aWStar )
        {
            MORIS_ERROR( false,
                    "IQI_Strong_Residual_Incompressible_NS::compute_dQIdu - not implemented\n." );
        }

        //------------------------------------------------------------------------------

        void
        IQI_SP_Crosswind_Incompressible_NS::compute_dQIdu(
                Vector< MSI::Dof_Type >& aDofType,
                Matrix< DDRMat >&        adQIdu )
        {
            MORIS_ERROR( false,
                    "IQI_Strong_Residual_Incompressible_NS::compute_dQIdu - not implemented\n." );
        }

        //------------------------------------------------------------------------------

        void
        IQI_SP_Crosswind_Incompressible_NS::compute_residual_strong_form_momentum( Matrix< DDRMat >& aRM )
        {
            // get the velocity and pressure FIs
            // FIXME protect dof type
            Field_Interpolator* tVelocityFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the properties
            const std::shared_ptr< Property >& tGravityProp = get_leader_property(IQI_Property_Type::GRAVITY);

            const std::shared_ptr< Property >& tThermalExpProp = get_leader_property(IQI_Property_Type::THERMAL_EXPANSION);

            const std::shared_ptr< Property >& tRefTempProp = get_leader_property(IQI_Property_Type::REF_TEMP);

            const std::shared_ptr< Property >& tInvPermeabProp = get_leader_property(IQI_Property_Type::INV_PERMEABILITY);

            // get the body load property
            const std::shared_ptr< Property >& tLoadProp = get_leader_property(IQI_Property_Type::BODY_LOAD);

            // get the incompressible fluid constitutive model
            const std::shared_ptr< Constitutive_Model >& tIncFluidCM = get_leader_constitutive_model(IQI_Constitutive_Type::INCOMPRESSIBLE_FLUID);

            // get the density property from CM
            const std::shared_ptr< Property >& tDensityProp = tIncFluidCM->get_property( "Density" );

            // get the density value
            real tDensity = tDensityProp->val()( 0 );

            // compute the residual strong form of momentum equation
            aRM = tDensity * trans( tVelocityFI->gradt( 1 ) )
                + tDensity * trans( tVelocityFI->gradx( 1 ) ) * tVelocityFI->val()
                - tIncFluidCM->divflux();

            // if body load
            if ( tLoadProp != nullptr )
            {
                // add contribution of body load term to momentum residual
                aRM -= tLoadProp->val();
            }

            // if permeability
            if ( tInvPermeabProp != nullptr )
            {
                // add Brinkman term to residual strong form
                aRM += tInvPermeabProp->val()( 0 ) * tVelocityFI->val();
            }

            // if gravity
            if ( tGravityProp != nullptr )
            {
                // add gravity to residual strong form
                aRM += tDensity * tGravityProp->val();

                // if thermal expansion and reference temperature
                if ( tThermalExpProp != nullptr && tRefTempProp != nullptr )
                {
                    // get the temperature field interpolator
                    // FIXME protect FI
                    Field_Interpolator* tTempFI = get_leader_fi_manager()->get_field_interpolators_for_type( MSI::Dof_Type::TEMP );

                    // add contribution to residual
                    aRM -= tDensity * tGravityProp->val() * tThermalExpProp->val()
                         * ( tTempFI->val() - tRefTempProp->val() );
                }
            }
        }

        //------------------------------------------------------------------------------
    }    // namespace fem
}    // namespace moris
