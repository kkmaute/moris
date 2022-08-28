/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.cpp
 *
 */

#include "cl_FEM_IWG_Diffusion_Virtual_Work_Ghost.hpp"
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

        IWG_Diffusion_Virtual_Work_Ghost::IWG_Diffusion_Virtual_Work_Ghost()
        {
            // set ghost flag
            mIsGhost = true;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize(  static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostVWOrder" ] = static_cast< uint >( IWG_Stabilization_Type::GHOST_VW );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif
            // set interpolation order
            IWG::set_interpolation_order();

            // get master index for residual dof type, indices for assembly
            uint tDofIndexMaster      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tDofIndexSlave      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tMasterFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tSlaveFI  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the diffusion constitutive model for master and slave
            const std::shared_ptr< Constitutive_Model > & tCMMasterDiff =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveDiff =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPGhost =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for master and slave
            real tMasterConductivity = tCMMasterDiff->constitutive()( 0, 0 );
            real tSlaveConductivity  = tCMSlaveDiff->constitutive()( 0, 0 );

            // get sub-matrices
            auto tResMaster = mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex } );

            auto tResSlave = mSet->get_residual()( 0 )(
                    { tSlaveResStartIndex, tSlaveResStopIndex } );

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, iOrder );

                // multiply common terms (note: needs to be stored as matrix
                Matrix< DDRMat > tPreMultiply =
                        tSPGhost->val()( 0 ) * trans( tNormalMatrix ) * tNormalMatrix *
                        ( tMasterConductivity * tMasterFI->gradx( iOrder ) - tSlaveConductivity * tSlaveFI->gradx( iOrder ) );

                // compute master residual
                tResMaster += aWStar * (
                        trans( tMasterFI->dnNdxn( iOrder ) ) * tPreMultiply );

                // compute slave residual
                tResSlave -= aWStar * (
                        trans( tSlaveFI->dnNdxn( iOrder ) ) * tPreMultiply );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Virtual_Work_Ghost::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // set interpolation order
            IWG::set_interpolation_order();

            // get master index for residual dof type, indices for assembly
            uint tDofIndexMaster      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tDofIndexSlave      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tMasterFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tSlaveFI  =
                    mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the diffusion constitutive model for master and slave
            const std::shared_ptr< Constitutive_Model > & tCMMasterDiff =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            const std::shared_ptr< Constitutive_Model > & tCMSlaveDiff =
                    mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSPGhost =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for master and slave
            real tMasterConductivity = tCMMasterDiff->constitutive()( 0, 0 );
            real tSlaveConductivity  = tCMSlaveDiff->constitutive()( 0, 0 );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, iOrder );

                // compute the Jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 );

                    auto tJacMaster = mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tDepStartIndex,       tDepStopIndex } );

                    auto tJacSlave =  mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tDepStartIndex,      tDepStopIndex } );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // add contribution to Jacobian
                        tJacMaster += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tMasterConductivity * tMasterFI->dnNdxn( iOrder ) );

                       tJacSlave -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tMasterConductivity * tMasterFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMMasterDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacMaster -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tMasterFI->gradx( iOrder ) * tCMMasterDiff->dConstdDOF( tDofType ) );

                        tJacSlave += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tMasterFI->gradx( iOrder ) * tCMMasterDiff->dConstdDOF( tDofType ) );
                    }
                }

                // compute the Jacobian for indirect dof dependencies through slave constitutive models
                uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    const Cell< MSI::Dof_Type > & tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                    // get index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 );

                    // get sub-matrices
                    auto tJacMaster = mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tDepStartIndex,       tDepStopIndex } );

                    auto tJacSlave = mSet->get_jacobian()(
                            { tSlaveResStartIndex, tSlaveResStopIndex },
                            { tDepStartIndex,      tDepStopIndex } );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // add contribution to Jacobian
                        tJacMaster -= aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tSlaveConductivity * tSlaveFI->dnNdxn( iOrder ) );

                        tJacSlave += aWStar * (
                                        tSPGhost->val()( 0 )
                                        * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                        * tNormalMatrix * tSlaveConductivity * tSlaveFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMSlaveDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJacMaster -= aWStar * (
                                tSPGhost->val()( 0 )
                                * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                * tNormalMatrix * tSlaveFI->gradx( iOrder ) * tCMSlaveDiff->dConstdDOF( tDofType ) );

                        tJacSlave += aWStar * (
                                tSPGhost->val()( 0 )
                                * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                * tNormalMatrix * tSlaveFI->gradx( iOrder ) * tCMSlaveDiff->dConstdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix(
                Matrix< DDRMat > & aFlatNormal,
                uint               aOrder )
        {
            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // switch on the ghost order
            switch( aOrder )
            {
                case 1 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        case 3 :
                        {
                            aFlatNormal = trans( mNormal );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 2 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 2, 3, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 2 ) = mNormal( 0 );

                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 1, 3 ) = mNormal( 2 );
                            aFlatNormal( 2, 3 ) = mNormal( 1 );

                            aFlatNormal( 0, 4 ) = mNormal( 2 );
                            aFlatNormal( 2, 4 ) = mNormal( 0 );

                            aFlatNormal( 0, 5 ) = mNormal( 1 );
                            aFlatNormal( 1, 5 ) = mNormal( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                case 3 :
                {
                    switch ( tSpaceDim )
                    {
                        case 2 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 3, 4, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );

                            aFlatNormal( 0, 2 ) = mNormal( 1 );
                            aFlatNormal( 1, 3 ) = mNormal( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 2, 2 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 2, 3 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        case 3 :
                        {
                            // set the normal matrix size
                            aFlatNormal.set_size( 6, 10, 0.0 );

                            aFlatNormal( 0, 0 ) = mNormal( 0 );
                            aFlatNormal( 1, 1 ) = mNormal( 1 );
                            aFlatNormal( 2, 2 ) = mNormal( 2 );

                            aFlatNormal( 0, 3 ) = mNormal( 1 );
                            aFlatNormal( 0, 4 ) = mNormal( 2 );

                            aFlatNormal( 1, 5 ) = mNormal( 0 );
                            aFlatNormal( 1, 6 ) = mNormal( 2 );

                            aFlatNormal( 2, 7 ) = mNormal( 0 );
                            aFlatNormal( 2, 8 ) = mNormal( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            aFlatNormal( 3, 3 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 3, 5 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 3, 9 ) = tSqrtOf2 * mNormal( 2 );

                            aFlatNormal( 4, 6 ) = tSqrtOf2 * mNormal( 1 );
                            aFlatNormal( 4, 8 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 4, 9 ) = tSqrtOf2 * mNormal( 0 );

                            aFlatNormal( 5, 4 ) = tSqrtOf2 * mNormal( 0 );
                            aFlatNormal( 5, 7 ) = tSqrtOf2 * mNormal( 2 );
                            aFlatNormal( 5, 9 ) = tSqrtOf2 * mNormal( 1 );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - Spatial dimensions can only be 2, 3." );
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Diffusion_Virtual_Work_Ghost::get_flat_normal_matrix - order not supported" );
                }
            }
        }
        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */

