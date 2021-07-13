/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.cpp
 *
 *  Created on: Feb 18, 2020
 *      Author: ritzert
 */
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Nitsche.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        IWG_Isotropic_Struc_Linear_Contact_Nitsche::IWG_Isotropic_Struc_Linear_Contact_Nitsche( sint aBeta )
        {
            // sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            // .resize: gives aValue:(The value to initialize the new elements with) and aCount:(new size of the Cell)
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "NitscheInterface" ] = static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_residual( real aWStar )
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
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            const std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tMasterWeight = tSPNitsche->val()( 1 );
            const real tSlaveWeight  = tSPNitsche->val()( 2 );

            // normal projection operator
            const  Matrix< DDRMat > tNormalProjector  = mNormal * trans( mNormal );

            // compute the jump
            const  Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute projection of displacement jump onto normal
            // const real tINormalJump = dot( tJump, mNormal );

            // evaluate average traction
            const Matrix< DDRMat > tTraction =
                    tMasterWeight * tCMMasterElasticity->traction( mNormal ) +
                    tSlaveWeight  * tCMSlaveElasticity->traction( mNormal );

            // compute contact pressure
            const real tIfcPressure = dot( tTraction, mNormal );

            // check for penetration (tINormalJump > 0 )
            // if ( tINormalJump >= 0.0 || tIfcPressure <= 0 )
            if ( tIfcPressure <= 0 )
            {
                // compute master residual
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex } ) += aWStar * (
                                - tFIMaster->N_trans() * tNormalProjector * tTraction
                                + mBeta * tMasterWeight * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump
                                + tNitsche * tFIMaster->N_trans() * tNormalProjector * tJump );

                // compute slave residual
                mSet->get_residual()( 0 )(
                        { tSlaveResStartIndex, tSlaveResStopIndex } ) += aWStar * (
                                + tFISlave->N_trans() * tNormalProjector * tTraction
                                + mBeta * tSlaveWeight * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump
                                - tNitsche * tFISlave->N_trans() * tNormalProjector * tJump );

                // check for nan, infinity
                MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                        "IWG_Isotropic_Struc_Linear_Interface::compute_residual - Residual contains NAN or INF, exiting!");
            }
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            const uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            const uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            const uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            const uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            const uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            const uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the elasticity constitutive model
            const std::shared_ptr< Constitutive_Model > & tCMMasterElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            const std::shared_ptr< Constitutive_Model > & tCMSlaveElasticity =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get the Nitsche stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            const real tNitsche      = tSPNitsche->val()( 0 );
            const real tMasterWeight = tSPNitsche->val()( 1 );
            const real tSlaveWeight  = tSPNitsche->val()( 2 );

            // normal projection operator
            const  Matrix< DDRMat > tNormalProjector  = mNormal * trans( mNormal );

            // compute the jump
            const  Matrix< DDRMat > tJump = tFIMaster->val() - tFISlave->val();

            // compute projection of displacement jump onto normal
            // const real tINormalJump = dot( tJump, mNormal );

            // evaluate average traction
            const Matrix< DDRMat > tTraction =
                    tMasterWeight * tCMMasterElasticity->traction( mNormal ) +
                    tSlaveWeight  * tCMSlaveElasticity->traction( mNormal );

            // compute contact pressure
            const real tIfcPressure = dot( tTraction, mNormal );

            // check for penetration (tINormalJump > 0 )
            // if ( tINormalJump >= 0.0 || tIfcPressure <= 0 )
            if ( tIfcPressure <= 0 )
            {
                // get number of master dof dependencies
                const uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

                // compute the Jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    const sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    const uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                    const uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                    // extract sub-matrices
                    auto tJacMM = mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } );

                    auto tJacSM =  mSet->get_jacobian()(
                            { tSlaveResStartIndex,  tSlaveResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } );

                    // compute Jacobian direct dependencies
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        tJacMM += aWStar * (
                                + mBeta * tMasterWeight * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tFIMaster->N()
                                + tNitsche * tFIMaster->N_trans() * tNormalProjector * tFIMaster->N() );

                        tJacSM += aWStar * (
                                + mBeta * tSlaveWeight * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tFIMaster->N()
                                - tNitsche * tFISlave->N_trans() * tNormalProjector * tFIMaster->N() );
                    }

                    // if dependency on the dof type
                    if ( tCMMasterElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacMM += aWStar * (
                                - tMasterWeight * tFIMaster->N_trans() * tNormalProjector * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal )
                                + mBeta * tMasterWeight * tCMMasterElasticity->dTestTractiondDOF( tDofType, mNormal, tNormalProjector * tJump, mResidualDofType( 0 ) ) );

                        tJacSM += aWStar * (
                                tMasterWeight * tFISlave->N_trans() * tNormalProjector * tCMMasterElasticity->dTractiondDOF( tDofType, mNormal ) );
                    }

                    // if dependency of stabilization parameters on the dof type
                    if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                    {
                        // get the derivatives of the SPs
                        const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 0 );
                        const Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 1 );
                        const Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdMasterDOF( tDofType ).get_row( 2 );

                        // get traction derivative
                        const Matrix< DDRMat > tTractionDer =
                                tCMMasterElasticity->traction( mNormal ) * tMasterWeightDer +
                                tCMSlaveElasticity->traction( mNormal ) * tSlaveWeightDer;

                        // add contribution to Jacobian
                        tJacMM += aWStar * (
                                - tFIMaster->N_trans() * tNormalProjector * tTractionDer +
                                + mBeta * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump * tMasterWeightDer
                                + tFIMaster->N_trans() * tNormalProjector * tJump * tNitscheDer );

                        tJacSM += aWStar * (
                                + tFISlave->N_trans() * tNormalProjector * tTractionDer
                                + mBeta * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump * tSlaveWeightDer
                                - tFISlave->N_trans() * tNormalProjector * tJump * tNitscheDer );
                    }
                }

                // compute the Jacobian for indirect dof dependencies through slave constitutive models
                uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    const Cell< MSI::Dof_Type > & tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    const sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                    const uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                    const uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                    // extract sub-matrices
                    auto tJacMS = mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tSlaveDepStartIndex,  tSlaveDepStopIndex  } );

                    auto tJacSS = mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tSlaveDepStartIndex, tSlaveDepStopIndex } );

                    // if dof type is residual dof type
                    if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        tJacMS += aWStar * (
                                - mBeta * tMasterWeight * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tFISlave->N()
                                - tNitsche * tFIMaster->N_trans() * tNormalProjector * tFISlave->N() );

                        tJacSS += aWStar * (
                                - mBeta * tSlaveWeight * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tFISlave->N()
                                + tNitsche * tFISlave->N_trans() * tNormalProjector * tFISlave->N() );
                    }

                    // if dependency on the dof type
                    if ( tCMSlaveElasticity->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacMS += aWStar * (
                                - tSlaveWeight * tFIMaster->N_trans() * tNormalProjector * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal ) );

                        tJacSS += aWStar * (
                                + tSlaveWeight * tFISlave->N_trans() * tNormalProjector * tCMSlaveElasticity->dTractiondDOF( tDofType, mNormal )
                                + mBeta * tSlaveWeight * tCMSlaveElasticity->dTestTractiondDOF( tDofType, mNormal, tNormalProjector * tJump, mResidualDofType( 0 ) ) );
                    }

                    // if dependency of stabilization parameters on the dof type
                    if ( tSPNitsche->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                    {
                        // get the derivatives of the SPs
                        const Matrix< DDRMat > tNitscheDer      = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 0 );
                        const Matrix< DDRMat > tMasterWeightDer = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 1 );
                        const Matrix< DDRMat > tSlaveWeightDer  = tSPNitsche->dSPdSlaveDOF( tDofType ).get_row( 2 );

                        // get traction derivative
                        const Matrix< DDRMat > tTractionDer =
                                tCMMasterElasticity->traction( mNormal ) * tMasterWeightDer +
                                tCMSlaveElasticity->traction( mNormal ) * tSlaveWeightDer;

                        // add contribution to Jacobian
                        tJacMS += aWStar * (
                                - tFIMaster->N_trans() * tNormalProjector * tTractionDer
                                + mBeta * tCMMasterElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump * tMasterWeightDer
                                + tFIMaster->N_trans() * tNormalProjector * tJump * tNitscheDer );

                        tJacSS += aWStar * (
                                tFISlave->N_trans() * tTractionDer
                                + mBeta * tCMSlaveElasticity->testTraction( mNormal, mResidualDofType( 0 ) ) * tNormalProjector * tJump * tSlaveWeightDer
                                - tFISlave->N_trans() * tNormalProjector * tJump * tNitscheDer );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Nitsche::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
