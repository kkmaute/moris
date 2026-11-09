/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strain_Energy.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Strain_Energy::IQI_Strain_Energy()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Bedding" ]   = static_cast< uint >( IQI_Property_Type::BEDDING );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IQI_Property_Type::THICKNESS );
        mPropertyMap[ "Select" ]    = static_cast< uint >( IQI_Property_Type::SELECT );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST );
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > &aQI )
    {
        // get select property
        const std::shared_ptr< Property > &tPropSelect =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        real tSelect = ( tPropSelect != nullptr ) ? tPropSelect->val()( 0 ) : 1.0;

        // check whether select is zero, if so return
        if ( std::abs( tSelect ) < MORIS_REAL_EPS )
        {
            aQI.set_size( 1, 1, 0.0 );
            return;
        }

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::BEDDING ) );

        // evaluate the QI
        aQI = 0.5 * tSelect * trans( tCMElasticity->flux() ) * tCMElasticity->strain();

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // check that field interpolator exists
            MORIS_ASSERT( tDisplacementFI != nullptr,
                    "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

            // compute body load contribution
            aQI += 0.5 * tSelect * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val();
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_QI( real aWStar )
    {
        // get select property
        const std::shared_ptr< Property > &tPropSelect =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        real tSelect = ( tPropSelect != nullptr ) ? tPropSelect->val()( 0 ) : 1.0;

        // check whether select is zero, if so return
        if ( std::abs( tSelect ) < MORIS_REAL_EPS )
        {
            return;
        }

        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::BEDDING ) );

        // get thickness property
        const std::shared_ptr< Property > &tPropThickness =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * tSelect * ( 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain() );

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // check that field interpolator exists
            MORIS_ASSERT( tDisplacementFI != nullptr,
                    "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

            // compute body load contribution
            mSet->get_QI()( tQIIndex ) += aWStar * tSelect *    //
                                          ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val() );
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_dQIdu( real aWStar )
    {
        // get select property
        const std::shared_ptr< Property > &tPropSelect =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        real tSelect = ( tPropSelect != nullptr ) ? tPropSelect->val()( 0 ) : 1.0;

        // check whether select is zero, if so return
        if ( std::abs( tSelect ) < MORIS_REAL_EPS )
        {
            return;
        }

        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::BEDDING ) );

        // get thickness property
        const std::shared_ptr< Property > &tPropThickness =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // if elasticity CM depends on dof type
            if ( tCMElasticity->check_dof_dependency( tDofType ) )
            {
                // compute dQIdu
                mSet->get_residual()( tQIIndex )(
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=                                   //
                        aWStar * 0.5 * tSelect * (                                                           //
                                trans( tCMElasticity->dFluxdDOF( tDofType ) ) * tCMElasticity->strain() +    //
                                trans( tCMElasticity->dStraindDOF( tDofType ) ) * tCMElasticity->flux() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator *tDisplacementFI =
                        mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // compute bedding contribution - displacements
                if ( tDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=    //
                            aWStar * tSelect * ( trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
                }

                if ( tPropBedding->check_dof_dependency( tDofType ) )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=    //
                            aWStar * tSelect * ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                }
            }

            // if select property depends on dof type
            if ( tPropSelect != nullptr )
            {
                if ( tPropSelect->check_dof_dependency( tDofType ) )
                {
                    MORIS_ERROR( false, "IQI_Strain_Energy::compute_dQIdu - select property cannot depend on dof type." );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_dQIdu(
            Vector< MSI::Dof_Type > &aDofType,
            Matrix< DDRMat >        &adQIdu )
    {
        // initialize derivative
        adQIdu.fill( 0.0 );

        // get select property
        const std::shared_ptr< Property > &tPropSelect =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::SELECT ) );

        real tSelect = ( tPropSelect != nullptr ) ? tPropSelect->val()( 0 ) : 1.0;

        // check whether select is zero, if so return
        if ( std::abs( tSelect ) < MORIS_REAL_EPS )
        {
            return;
        }

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IQI_Property_Type::BEDDING ) );

        // if elasticity CM depends on dof type
        if ( tCMElasticity->check_dof_dependency( aDofType ) )
        {
            // compute dQIdu
            adQIdu = 0.5 * tSelect * (                                                                    //
                             trans( tCMElasticity->dFluxdDOF( aDofType ) ) * tCMElasticity->strain() +    //
                             trans( tCMElasticity->dStraindDOF( aDofType ) ) * tCMElasticity->flux() );
        }

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // compute bedding contribution - displacements
            if ( aDofType( 0 ) == MSI::Dof_Type::UX )
            {
                adQIdu += tSelect * trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val();
            }

            if ( tPropBedding->check_dof_dependency( aDofType ) )
            {
                adQIdu += 0.5 * tSelect * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( aDofType );
            }
        }

        // if select property depends on dof type
        if ( tPropSelect != nullptr )
        {
            if ( tPropSelect->check_dof_dependency( aDofType ) )
            {
                MORIS_ERROR( false, "IQI_Strain_Energy::compute_dQIdu - select property cannot depend on dof type." );
            }
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
