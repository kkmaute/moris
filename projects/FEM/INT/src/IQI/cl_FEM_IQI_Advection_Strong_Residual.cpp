/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Advection_Strong_Residual.cpp
 *
 */

#include "cl_FEM_IQI_Advection_Strong_Residual.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Advection_Strong_Residual::IQI_Advection_Strong_Residual()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::ADVECTION_STRONG_RESIDUAL;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Load" ] = static_cast< uint >( IQI_Property_Type::BODY_LOAD );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IQI_Constitutive_Type::DIFFUSION );
    }

    //------------------------------------------------------------------------------

    void IQI_Advection_Strong_Residual::compute_QI( Matrix< DDRMat >& aQI )
    {
        Field_Interpolator* tFIVelocity =
                mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get the diffusion CM
        const std::shared_ptr< Constitutive_Model >& tCMDiffusion =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::DIFFUSION ) );

        // get body load property
        const std::shared_ptr< Property >& tPropLoad =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::BODY_LOAD ) );

        Matrix< DDRMat > tRT = tCMDiffusion->EnergyDot() + tFIVelocity->val_trans() * tCMDiffusion->gradEnergy() - tCMDiffusion->divflux();

        // if body load exists
        if ( tPropLoad != nullptr )
        {
            // add energy source term
            tRT -= tPropLoad->val();
        }

        // evaluate the QI
        aQI = tRT( 0, 0 );
    }

    //------------------------------------------------------------------------------

    void IQI_Advection_Strong_Residual::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate strong form
        Matrix< DDRMat > tQI( 1, 1 );
        this->compute_QI( tQI );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * tQI;
    }

    //------------------------------------------------------------------------------

    void IQI_Advection_Strong_Residual::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false,
                "IQI_Advection_Strong_Residual::compute_dQIdu - not implemented\n." );
    }

    //------------------------------------------------------------------------------

    void IQI_Advection_Strong_Residual::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false,
                "IQI_Advection_Strong_Residual::compute_dQIdu - not implemented\n." );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
