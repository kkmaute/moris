
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        IWG_Isotropic_Struc_Linear_Bulk::IWG_Isotropic_Struc_Linear_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ] = IWG_Property_Type::LOAD;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::set_property(
                std::shared_ptr< Property > aProperty,
                std::string                 aPropertyString,
                mtk::Master_Slave           aIsMaster )
        {
            // check that aPropertyString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Isotropic_Struc_Linear_Bulk::set_property - Unknown aPropertyString: " ) +
                    aPropertyString;
            MORIS_ERROR( mPropertyMap.find( aPropertyString ) != mPropertyMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Isotropic_Struc_Linear_Bulk::set_property - No slave allowed." );

            // set the property in the property cell
            this->get_properties( aIsMaster )( static_cast< uint >( mPropertyMap[ aPropertyString ] ) ) = aProperty;
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster  )
        {
            // check that aConstitutiveString makes sense
            std::string tErrMsg =
                    std::string( "IWG_Isotropic_Struc_Linear_Bulk::set_constitutive_model - Unknown aConstitutiveString: " ) +
                    aConstitutiveString;
            MORIS_ERROR( mConstitutiveMap.find( aConstitutiveString ) != mConstitutiveMap.end(), tErrMsg.c_str() );

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Isotropic_Struc_Linear_Bulk::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get body load property
            std::shared_ptr< Property > tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get elasticity CM
            std::shared_ptr< Constitutive_Model > tCMElasticty =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // compute the residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                    += aWStar * ( trans( tCMElasticty->testStrain() ) * tCMElasticty->flux() );

            // if body load
            if ( tPropLoad != nullptr )
            {
                // compute body load contribution
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                        -= aWStar * ( trans( tDisplacementFI->N() ) * tPropLoad->val() );
            }

        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            Field_Interpolator * tDisplacementFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ));

            // get body load property
            std::shared_ptr< Property > tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get elasticity CM
            std::shared_ptr< Constitutive_Model > tCMElasticty =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the number of master dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the master dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if body load
                if ( tPropLoad != nullptr )
                {
                    // if property depends on the dof type
                    if ( tPropLoad->check_dof_dependency( tDofType ) )
                    {
                        // compute the contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tMasterDepStartIndex, tMasterDepStopIndex } ) -= aWStar * (
                                        trans( tDisplacementFI->N() ) * tPropLoad->dPropdDOF( tDofType ) );
                    }
                }

                // if constitutive model depends on the dof type
                if ( tCMElasticty->check_dof_dependency( tDofType ) )
                {
                    // compute the contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tCMElasticty->testStrain() ) * tCMElasticty->dFluxdDOF( tDofType ) );
                }
            }
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
