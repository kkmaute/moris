/*
 * cl_FEM_IWG_Compressible_NS_Density_Bulk.cpp
 *
 *  Created on: Jul 28, 2020
 *      Author: wunsch
 */

#include "cl_FEM_IWG_Compressible_NS_Density_Bulk.hpp"

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------
        IWG_Compressible_NS_Density_Bulk::IWG_Compressible_NS_Density_Bulk()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Fluid" ] = IWG_Constitutive_Type::FLUID;
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Compressible_NS_Density_Bulk::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Compressible_NS_Density_Bulk::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster  )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Compressible_NS_Density_Bulk::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Compressible_NS_Density_Bulk::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::set_stabilization_parameter - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tDensityFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX ); //FIXME this need to be protected

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tFluidCM =  mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // compute the residual weak form
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } ) += aWStar * (
                                trans( tDensityFI->N() ) * tDensityFI->gradt( 1 )
                                - trans( tDensityFI->dnNdxn( 1 ) ) * tDensityFI->val()( 0 ) * tVelocityFI->val() );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_jacobian( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here pressure), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the field interpolators
            Field_Interpolator * tDensityFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );
            Field_Interpolator * tVelocityFI = mMasterFIManager->get_field_interpolators_for_type( mDofVelocity );

            // get the constitutive model
            std::shared_ptr< Constitutive_Model > tFluidCM =  mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID ) );

            // compute the jacobian for dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is density, add diagonal term (density - density DoF types)
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // add contibution
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tDensityFI->N() ) * tDensityFI->dnNdxn( 1 )  -
                                    trans( tDensityFI->dnNdxn( 1 ) ) * tVelocityFI->val() * tDensityFI->N() );
                }

                // if dof type is velocity, add the mixed term (velocity - density DoF types)
                if( tDofType( 0 ) == mDofVelocity )
                {
                    // compute the jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - 1.0 * trans( tDensityFI->dnNdxn( 1 ) ) * tDensityFI->val()( 0 ) * tVelocityFI->N() );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_residual_strong_form(
                Matrix< DDRMat > & aRM,
                real             & aRC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_residual_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void IWG_Compressible_NS_Density_Bulk::compute_jacobian_strong_form(
                moris::Cell< MSI::Dof_Type >   aDofTypes,
                Matrix< DDRMat >             & aJM,
                Matrix< DDRMat >             & aJC )
        {
            MORIS_ERROR( false, "IWG_Compressible_NS_Density_Bulk::compute_jacobian_strong_form - Not implemented." );
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
