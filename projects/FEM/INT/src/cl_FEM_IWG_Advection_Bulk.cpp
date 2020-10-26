
#include "cl_FEM_IWG_Advection_Bulk.hpp"
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

        IWG_Advection_Bulk::IWG_Advection_Bulk()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = IWG_Constitutive_Type::DIFFUSION;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = IWG_Stabilization_Type::SUPG;
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::set_constitutive_model(
                std::shared_ptr< Constitutive_Model > aConstitutiveModel,
                std::string                           aConstitutiveString,
                mtk::Master_Slave                     aIsMaster )
        {
            // check that aConstitutiveString makes sense
            if ( mConstitutiveMap.find( aConstitutiveString ) == mConstitutiveMap.end() )
            {
                std::string tErrMsg =
                        std::string("IWG_Advection_Bulk::set_constitutive_model - Unknown aConstitutiveString: ") +
                        aConstitutiveString;

                MORIS_ERROR( false, tErrMsg.c_str() );
            }

            // check no slave allowed
            MORIS_ERROR( aIsMaster == mtk::Master_Slave::MASTER,
                    "IWG_Advection_Bulk::set_constitutive_model - No slave allowed." );

            // set the constitutive model in the constitutive model cell
            this->get_constitutive_models( aIsMaster )( static_cast< uint >( mConstitutiveMap[ aConstitutiveString ] ) ) = aConstitutiveModel;
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::set_stabilization_parameter(
                std::shared_ptr< Stabilization_Parameter > aStabilizationParameter,
                std::string                                aStabilizationString )
        {
            // check that aStabilizationString makes sense
            std::string tErrMsg =
                    std::string( "IWG_GGLS_Diffusion_Phase_Change_Bulk::set_stabilization_parameter - Unknown aStabilizationString: " ) +
                    aStabilizationString;
            MORIS_ERROR( mStabilizationMap.find( aStabilizationString ) != mStabilizationMap.end(), tErrMsg.c_str() );


            // set the stabilization parameter in the stabilization parameter cell
            this->get_stabilization_parameters()( static_cast< uint >( mStabilizationMap[ aStabilizationString ] ) ) = aStabilizationParameter;
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute the residual strong form
            Matrix< DDRMat > tRT;
            this->compute_residual_strong_form( tRT );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            trans( tFITemp->N() ) * trans( tFIVelocity->val() ) * tCMDiffusion->gradEnergy() +
                            tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tRT( 0, 0 ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Advection_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute the residual strong form
            Matrix< DDRMat > tRT;
            this->compute_residual_strong_form( tRT );

            // get number of master dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if velocity dof type
                // FIXME protect dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFITemp->N() ) * trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N() +
                                    tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->N() * tRT( 0, 0 ) );
                }

                // if diffusion CM depends on dof type
                if( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar *
                            ( trans( tFITemp->N() ) * trans( tFIVelocity->val() ) * tCMDiffusion->dGradEnergydDOF( tDofType ) );
                }

                // compute the jacobian strong form
                Matrix< DDRMat > tJT;
                this->compute_jacobian_strong_form( tDofType, tJT );

                // compute the jacobian contribution from strong form
                mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tJT );

                // if SUPG stabilization parameter depends on dof type
                if( tSPSUPG->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tRT( 0, 0 ) * tSPSUPG->dSPdMasterDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Advection_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Advection_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Advection_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_residual_strong_form( Matrix< DDRMat > & aRT )
        {
            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            aRT = tCMDiffusion->EnergyDot() +
                    trans( tFIVelocity->val() ) * tCMDiffusion->gradEnergy() -
                    tCMDiffusion->divflux();
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian_strong_form (
                moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >             & aJT )
        {
            // get the res dof and the derivative dof FIs
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // initialize aJT
            aJT.set_size( 1, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

            // get the diffusion CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // if CM diffusion depends on dof type
            if( tCMDiffusion->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to jacobian strong form
                aJT +=
                        tCMDiffusion->dEnergyDotdDOF( aDofTypes ) +
                        trans( tFIVelocity->val() ) * tCMDiffusion->dGradEnergydDOF( aDofTypes ) -
                        tCMDiffusion->ddivfluxdu( aDofTypes );
            }

            // if derivative wrt to velocity dof type
            if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                aJT +=
                        trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N() ;
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
