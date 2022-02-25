// FEM/INT/src
#include "cl_FEM_IWG_Spalart_Allmaras_Turbulence_Bulk.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

        //------------------------------------------------------------------------------

        IWG_Spalart_Allmaras_Turbulence_Bulk::IWG_Spalart_Allmaras_Turbulence_Bulk()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "SpalartAllmarasTurbulence" ] = static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "SUPG" ] = static_cast< uint >( IWG_Stabilization_Type::SUPG );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual( real aWStar )
        {
            // check master field interpolators
#ifdef DEBUG
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity = tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute residual of the strong form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get sub-matrix of residual
            auto tRes = mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex } );

            // compute the residual weak form
            tRes += aWStar
                  * ( tFIViscosity->N_trans()
                                  * ( tFIViscosity->gradt( 1 ) + trans( tModVelocity ) * tFIViscosity->gradx( 1 )
                                          - tCMSATurbulence->production_term() + tCMSATurbulence->wall_destruction_term() )
                          + trans( tFIViscosity->dnNdxn( 1 ) )
                                    * ( tCMSATurbulence->diffusion_coefficient()( 0 ) * tFIViscosity->gradx( 1 )
                                            + tModVelocity * tSPSUPG->val()( 0 ) * tR( 0 ) ) );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here velocity), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the SUPG stabilization parameter
            const std::shared_ptr< Stabilization_Parameter >& tSPSUPG =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::SUPG ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity = tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute residual of the strng form
            Matrix< DDRMat > tR;
            this->compute_residual_strong_form( tR );

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the dof dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Cell< MSI::Dof_Type >& tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                const uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                const uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // extract sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tMasterResStartIndex, tMasterResStopIndex }, { tMasterDepStartIndex, tMasterDepStopIndex } );

                // if residual dof type (here viscosity)
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute dModVelocitydModViscosity
                    Matrix< DDRMat > tModVelocityDer = -mCb2 * tFIViscosity->dnNdxn( 1 ) / mSigma;

                    // compute the jacobian
                    tJac += aWStar
                          * ( tFIViscosity->N_trans()
                                          * ( tFIViscosity->dnNdtn( 1 )
                                                  + trans( tFIVelocity->val() - 2.0 * mCb2 * tFIViscosity->gradx( 1 ) / mSigma )
                                                            * tFIViscosity->dnNdxn( 1 ) )
                                  + trans( tFIViscosity->dnNdxn( 1 ) )
                                            * ( tCMSATurbulence->diffusion_coefficient()( 0 ) * tFIViscosity->dnNdxn( 1 )
                                                    + tModVelocityDer * tSPSUPG->val()( 0 ) * tR( 0 ) ) );
                }
                
                // if velocity dof type
                // FIXME protect dof type
                else if ( tDofType( 0 ) == MSI::Dof_Type::VX )
                {
                    // compute dModVelocitydVelocity
                    const Matrix< DDRMat >& tModVelocityDer = tFIVelocity->N();

                    // compute the jacobian
                    tJac += aWStar
                          * ( tFIViscosity->N_trans() * trans( tFIViscosity->gradx( 1 ) ) * tModVelocityDer
                                  + trans( tFIViscosity->dnNdxn( 1 ) ) * tModVelocityDer * tSPSUPG->val()( 0 ) * tR( 0 ) );
                }

                // if turbulence model depends on dof
                if ( tCMSATurbulence->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    tJac += aWStar
                          * ( tFIViscosity->N_trans()
                                          * ( -1.0 * tCMSATurbulence->dproductiontermdu( tDofType )
                                                  + tCMSATurbulence->dwalldestructiontermdu( tDofType ) )
                                  + trans( tFIViscosity->dnNdxn( 1 ) ) * tFIViscosity->gradx( 1 )
                                            * tCMSATurbulence->ddiffusioncoeffdu( tDofType ) );
                }

                // if SP SUPG depends on dof
                if ( tSPSUPG->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * ( trans( tFIViscosity->dnNdxn( 1 ) ) * tModVelocity * tR( 0 ) * tSPSUPG->dSPdMasterDOF( tDofType ) );
                }

                // compute jacobian of the strong form
                Matrix< DDRMat > tJ;
                this->compute_jacobian_strong_form( tDofType, tJ );

                // compute the jacobian
                tJac += aWStar * ( trans( tFIViscosity->dnNdxn( 1 ) ) * tModVelocity * tSPSUPG->val()( 0 ) * tJ );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_and_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------
        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_dRdp( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_residual_strong_form( Matrix< DDRMat >& aR )
        {
            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // compute modified velocity
            Matrix< DDRMat > tModVelocity = tFIVelocity->val() - mCb2 * tFIViscosity->gradx( 1 ) / mSigma;

            // compute strong form of residual
            aR = tFIViscosity->gradt( 1 ) + trans( tModVelocity ) * tFIViscosity->gradx( 1 ) - tCMSATurbulence->production_term()
               + tCMSATurbulence->wall_destruction_term() - tCMSATurbulence->divflux();
        }

        //------------------------------------------------------------------------------

        void
        IWG_Spalart_Allmaras_Turbulence_Bulk::compute_jacobian_strong_form(
                const moris::Cell< MSI::Dof_Type >& aDofTypes,
                Matrix< DDRMat >&                   aJ )
        {
            // get the residual viscosity FI
            Field_Interpolator* tFIViscosity = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the velocity FI
            // FIXME protect dof type
            Field_Interpolator* tFIVelocity = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

            // get the SA turbulence CM
            const std::shared_ptr< Constitutive_Model >& tCMSATurbulence =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::SPALART_ALLMARAS_TURBULENCE ) );

            // get the der FI
            Field_Interpolator* tFIDer = mMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init aJ
            aJ.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // if dof type is residual dof type (here viscosity)
            if ( aDofTypes( 0 ) == mResidualDofType( 0 )( 0 ) )
            {
                aJ = tFIViscosity->dnNdtn( 1 )
                   + trans( tFIVelocity->val() - 2.0 * mCb2 * tFIViscosity->gradx( 1 ) / mSigma ) * tFIViscosity->dnNdxn( 1 );    // -
            }
            // if dof type is velocity dof type
            // FIXME protect dof type
            else if ( aDofTypes( 0 ) == MSI::Dof_Type::VX )
            {
                aJ = trans( tFIViscosity->gradx( 1 ) ) * tFIVelocity->N();
            }
            else
            {
                aJ.fill( 0.0 );
            }

            // if turbulence model depends on dof
            if ( tCMSATurbulence->check_dof_dependency( aDofTypes ) )
            {
                // add contribution to jacobian
                aJ += -1.0 * tCMSATurbulence->dproductiontermdu( aDofTypes )
                    + tCMSATurbulence->dwalldestructiontermdu( aDofTypes )
                    - tCMSATurbulence->ddivfluxdu( aDofTypes );
            }
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
