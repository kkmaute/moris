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

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IQI_Strain_Energy::IQI_Strain_Energy()
        {
            // set fem IQI type
            mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Bedding" ]   = static_cast< uint >( IWG_Property_Type::BEDDING );
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST );
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > & aQI )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // evaluate the QI
            aQI = 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain();

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator * tDisplacementFI =
                        mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // check that field interpolator exists
                MORIS_ASSERT( tDisplacementFI != nullptr,
                        "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist.");

                // compute body load contribution
                aQI += 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val();
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_QI( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get index for QI
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // evaluate the QI
            mSet->get_QI()( tQIIndex ) += aWStar * (
                    0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain() );

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator * tDisplacementFI =
                        mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // check that field interpolator exists
                MORIS_ASSERT( tDisplacementFI != nullptr,
                        "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist.");

                // compute body load contribution
                mSet->get_QI()( tQIIndex ) += aWStar * (
                        0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val() );
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu( real aWStar )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the column index to assemble in residual
            sint tQIIndex = mSet->get_QI_assembly_index( mName );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute dQIdu for indirect dof dependencies
            for( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDof );

                // get master index for residual dof type, indices for assembly
                uint tMasterDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
                uint tMasterDepStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

                // if elasticity CM depends on dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute dQIdu
                    mSet->get_residual()( tQIIndex )(
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * 0.5 * (
                                    trans( tCMElasticity->dFluxdDOF( tDofType ) )   * tCMElasticity->strain() +
                                    trans( tCMElasticity->dStraindDOF( tDofType ) ) * tCMElasticity->flux() );
                }

                // if bedding
                if ( tPropBedding != nullptr )
                {
                    // get field interpolator for displacements
                    Field_Interpolator * tDisplacementFI =
                            mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                    // compute bedding contribution - displacements
                    if( tDofType( 0 ) == MSI::Dof_Type::UX )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
                    }

                    if ( tPropBedding->check_dof_dependency( tDofType ) )
                    {
                        mSet->get_residual()( tQIIndex )(
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                    }
                }
            }
        }

        //------------------------------------------------------------------------------

        void IQI_Strain_Energy::compute_dQIdu(
                moris::Cell< MSI::Dof_Type > & aDofType,
                Matrix< DDRMat >             & adQIdu )
        {
            // check the point is inside the bounded box
            if ( !this->is_within_box_bounds() )
            {
                return;
            }

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // initialize derivative
            adQIdu.fill( 0.0 );

            // if elasticity CM depends on dof type
            if ( tCMElasticity->check_dof_dependency( aDofType ) )
            {
                // compute dQIdu
                adQIdu = 0.5 * (
                        trans( tCMElasticity->dFluxdDOF( aDofType ) )   * tCMElasticity->strain() +
                        trans( tCMElasticity->dStraindDOF( aDofType ) ) * tCMElasticity->flux() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator * tDisplacementFI =
                        mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // compute bedding contribution - displacements
                if( aDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    adQIdu += trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val();
                }

                if ( tPropBedding->check_dof_dependency( aDofType ) )
                {
                    adQIdu += 0.5 *
                            trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( aDofType ) ;
                }
            }
        }

        //------------------------------------------------------------------------------

        bool IQI_Strain_Energy::is_within_box_bounds()
        {
                // check if the box bounds are empty then skip
                if ( mParameters.empty() ) 
                {
                        return true;
                }

                //if the box bounds are not empty then check if it is inside the box
                else
                {
                        // get the coordinate 
                        const Matrix<DDRMat> & tGaussPoint = mMasterFIManager->get_IG_geometry_interpolator()->valx();

                        //check if the calculation point coordinates are more then lower corner of the box 
                        bool tLowerBound = std::equal(mParameters(0).begin(), mParameters(0).end(), tGaussPoint.begin(),tGaussPoint.end() , 
                        [](real aA, real aB) -> bool { return aA < aB ;} );
                        
                        //check if the calculation point coordinates are less then upper corner of the box 
                        bool tUpperBound = std::equal(tGaussPoint.begin(),tGaussPoint.end(), mParameters(1).begin(), mParameters(1).end(),
                         [](real aA, real aB) -> bool { return aA < aB ;} );
                        
                        //combine the two bounds that satisfy both
                        return tUpperBound and tLowerBound; 
                }
        }

        //------------------------------------------------------------------------------
    }/* end_namespace_fem */
}/* end_namespace_moris */

