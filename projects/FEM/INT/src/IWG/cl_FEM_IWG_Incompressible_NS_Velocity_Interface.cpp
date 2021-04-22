
#include "cl_FEM_IWG_Incompressible_NS_Velocity_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Incompressible_NS_Velocity_Interface::IWG_Incompressible_NS_Velocity_Interface( sint aBeta )
        {
            // set mBeta for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "IncompressibleFluid" ] = static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tMasterWeight * tCMMasterFluid->traction( mNormal ) +
                    tSlaveWeight  * tCMSlaveFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFIMaster->val() - tFISlave->val();

            // compute master residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - tFIMaster->N_trans() * tTractionFluid
                            - mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity
                            + tNitsche * tFIMaster->N_trans() * tJumpVelocity ) ;

            // compute slave residual
            mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            + tFISlave->N_trans() * tTractionFluid
                            - mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity
                            - tNitsche * tFISlave->N_trans() * tJumpVelocity );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Incompressible_NS_Velocity_Interface::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the fluid constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterFluid =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveFluid =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::FLUID_INCOMPRESSIBLE ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );
            real tNitsche      = tSPNitsche->val()( 0 );
            real tMasterWeight = tSPNitsche->val()( 1 );
            real tSlaveWeight  = tSPNitsche->val()( 2 );

            // evaluate average traction
            Matrix< DDRMat > tTractionFluid =
                    tMasterWeight * tCMMasterFluid->traction( mNormal ) +
                    tSlaveWeight  * tCMSlaveFluid->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJumpVelocity = tFIMaster->val() - tFISlave->val();

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIMaster->N()
                                    + tNitsche * tFIMaster->N_trans() * tFIMaster->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFIMaster->N()
                                    - tNitsche * tFISlave->N_trans() * tFIMaster->N() );
                }

                // if dependency on the dof type
                if ( tCMMasterFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tFIMaster->N_trans() * tMasterWeight * tCMMasterFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tMasterWeight * tCMMasterFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tJumpVelocity, mResidualDofType( 0 ) ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    + tFISlave->N_trans() * tMasterWeight * tCMMasterFluid->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            tCMMasterFluid->traction( mNormal ) * tMasterWeightDer +
                            tCMSlaveFluid->traction( mNormal )  * tSlaveWeightDer;

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - tFIMaster->N_trans() * tTractionDer
                                    - mBeta * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tMasterWeightDer
                                    + tFIMaster->N_trans() * tJumpVelocity * tNitscheDer );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    + tFISlave->N_trans() * tTractionDer
                                    - mBeta * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tSlaveWeightDer
                                    - tFISlave->N_trans() * tJumpVelocity * tNitscheDer );
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    + mBeta * tMasterWeight * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFISlave->N()
                                    - tNitsche * tFIMaster->N_trans() * tFISlave->N() );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    + mBeta * tSlaveWeight * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tFISlave->N()
                                    + tNitsche * tFISlave->N_trans() * tFISlave->N() );
                }

                // if dependency on the dof type
                if ( tCMSlaveFluid->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - tFIMaster->N_trans() * tSlaveWeight * tCMSlaveFluid->dTractiondDOF( tDofType, mNormal ) );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    + tFISlave->N_trans() * tSlaveWeight * tCMSlaveFluid->dTractiondDOF( tDofType, mNormal )
                                    - mBeta * tSlaveWeight * tCMSlaveFluid->dTestTractiondDOF(
                                            tDofType, mNormal, tJumpVelocity, mResidualDofType( 0 ) ) );
                }

                // if dependency of stabilization parameters on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // get the derivatives of the SPs
                    Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 0 );
                    Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                    Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                    // get traction derivative
                    Matrix< DDRMat > tTractionDer =
                            tCMMasterFluid->traction( mNormal ) * tMasterWeightDer +
                            tCMSlaveFluid->traction( mNormal )  * tSlaveWeightDer;

                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } ) += aWStar * (
                                    - tFIMaster->N_trans() * tTractionDer
                                    - mBeta * trans( tCMMasterFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tMasterWeightDer
                                    + tFIMaster->N_trans() * tJumpVelocity * tNitscheDer );

                    mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } ) += aWStar * (
                                    + tFISlave->N_trans() * tTractionDer
                                    - mBeta * trans( tCMSlaveFluid->testTraction( mNormal, mResidualDofType( 0 ) ) ) * tJumpVelocity * tSlaveWeightDer
                                    - tFISlave->N_trans() * tJumpVelocity * tNitscheDer );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Incompressible_NS_Velocity_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Incompressible_NS_Velocity_Interface::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Incompressible_NS_Velocity_Interface::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
