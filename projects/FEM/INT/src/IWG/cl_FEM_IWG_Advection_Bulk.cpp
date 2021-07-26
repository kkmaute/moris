
#include "cl_FEM_IWG_Advection_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Advection_Bulk::IWG_Advection_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ] = static_cast< uint >( IWG_Property_Type::BODY_LOAD );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFFUSION );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = static_cast< uint >( IWG_Stabilization_Type::SUPG );
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            const uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            const uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            const uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute the residual strong form
            Matrix< DDRMat > tRT;
            this->compute_residual_strong_form( tRT );

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } ) += aWStar * (
                            tFITemp->N_trans() * tFIVelocity->val_trans() * tCMDiffusion->gradEnergy() +
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
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual dof (here temperature) field interpolator
            Field_Interpolator* tFITemp = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute the residual strong form
            Matrix< DDRMat > tRT;
            this->compute_residual_strong_form( tRT );

            // get number of master dof dependencies
            const  uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                const uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                const uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto jac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if velocity dof type
                // FIXME protect dof type
                if( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // add contribution to Jacobian
                    jac += aWStar * (
                            tFITemp->N_trans() * trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N() +
                            tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->N() * tRT( 0, 0 ) );
                }

                // if diffusion CM depends on dof type
                if( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    jac += aWStar *
                            ( tFITemp->N_trans() * tFIVelocity->val_trans() * tCMDiffusion->dGradEnergydDOF( tDofType ) );
                }

                // compute the Jacobian strong form
                Matrix< DDRMat > tJT;
                this->compute_jacobian_strong_form( tDofType, tJT );

                // compute the Jacobian contribution from strong form
                jac += aWStar * (
                        tSPSUPG->val()( 0 ) * trans( tFITemp->dnNdxn( 1 ) ) * tFIVelocity->val() * tJT );

                // if SUPG stabilization parameter depends on dof type
                if( tSPSUPG->check_dof_dependency( tDofType ) )
                {
                    // add contribution to Jacobian
                    jac += aWStar * (
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
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            aRT = tCMDiffusion->EnergyDot() +
                    tFIVelocity->val_trans() * tCMDiffusion->gradEnergy() -
                    tCMDiffusion->divflux();

            // if body load exists
            if ( tPropLoad != nullptr )
            {
                // add energy source term
                aRT -= tPropLoad->val();
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Advection_Bulk::compute_jacobian_strong_form (
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & aJT )
        {
            // get the res dof and the derivative dof FIs
            Field_Interpolator * tFIDer =
                    mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity dof field interpolator
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity =
                    mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // initialize aJT
            aJT.set_size( 1, tFIDer->get_number_of_space_time_coefficients());

            // get the diffusion CM
            const std::shared_ptr< Constitutive_Model > & tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFFUSION ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::BODY_LOAD ) );

            // if CM diffusion depends on dof type
            if( tCMDiffusion->check_dof_dependency( aDofTypes ) )
            {
                // compute contribution to Jacobian strong form
                aJT =   tCMDiffusion->dEnergyDotdDOF( aDofTypes ) +
                        tFIVelocity->val_trans() * tCMDiffusion->dGradEnergydDOF( aDofTypes ) -
                        tCMDiffusion->ddivfluxdu( aDofTypes );
            }
            else
            {
                aJT.fill(0.0);
            }

            // if derivative wrt to velocity dof type
            if( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                aJT += trans( tCMDiffusion->gradEnergy() ) * tFIVelocity->N() ;
            }

            // if body load exists and depends on DOFs
            if ( tPropLoad != nullptr )
            {
                if ( tPropLoad->check_dof_dependency( aDofTypes ) )
                {
                    aJT -= tPropLoad->dPropdDOF( aDofTypes );
                }
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
