/*
 * cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.cpp
 *
 *  Created on: Jan 27, 2020
 *      Author: ritzert
 */
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Contact_Penalty.hpp"
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
        IWG_Isotropic_Struc_Linear_Contact_Penalty::IWG_Isotropic_Struc_Linear_Contact_Penalty()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "PenaltyContact" ] = IWG_Stabilization_Type::NITSCHE_INTERFACE;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // FIXME what is this parameter???
            moris::real tTempParam = 1e13;

            // evaluate the gap
            Matrix< DDRMat > tGap = trans( mNormal ) * ( tFISlave->val() - tFIMaster->val() ); // mNormal is normal on Master side
            if ( tGap( 0, 0 ) > 0 )
            {
                tGap( 0, 0 ) = 0;
            }

            // compute contact residual on master side
            mSet->get_residual()( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            -= aWStar * ( tTempParam * trans( tFIMaster -> N() ) * mNormal * tGap( 0, 0 ) );

            // compute contact residual on slave side
            mSet->get_residual()( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
            += aWStar * ( tTempParam * trans( tFISlave -> N() ) * mNormal * tGap( 0, 0 ) );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tSlaveDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tSlaveDofIndex )( 0, 1 );

            // get the master field interpolator for the residual dof type
            Field_Interpolator * tFIMaster = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the slave field interpolator for the residual dof type
            Field_Interpolator * tFISlave  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPPenalty
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::NITSCHE_INTERFACE ) );

            // create quadrature of normal
            Matrix< DDRMat > tNormQ = mNormal * trans ( mNormal );

            // FIXME what is this parameter???
            moris::real tTempParam = 1e13;

            // evaluate the gap
            Matrix< DDRMat > tGap = trans( mNormal ) * ( tFISlave->val() - tFIMaster->val() );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the master dof dependencies
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // compute jacobian direct dependencies
                if ( tDofType( 0 ) == mResidualDofType( 0 ) && tGap( 0, 0 ) <= 0 )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += aWStar * ( tTempParam * trans( tFIMaster -> N() ) * tNormQ * tFIMaster -> N() );

                    mSet->get_jacobian()( { tSlaveResStartIndex,  tSlaveResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    -= aWStar * ( tTempParam * trans( tFISlave -> N() ) * tNormQ * tFIMaster -> N() );
                }
            }

            // get number of slave dof dependencies
            uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();

            // loop over slave dof dependencies
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex        = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                uint tSlaveDepStartIndex = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 0 );
                uint tSlaveDepStopIndex  = mSet->get_jac_dof_assembly_map()( tSlaveDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 ) && tGap( 0, 0 ) <= 0 )
                {
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tSlaveDepStartIndex,  tSlaveDepStopIndex  } )
                    -= aWStar * ( tTempParam * trans( tFIMaster -> N() ) * tNormQ * tFISlave -> N() );

                    mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                          { tSlaveDepStartIndex, tSlaveDepStopIndex } )
                    += aWStar * ( tTempParam * trans( tFISlave -> N() ) * tNormQ * tFISlave -> N() );
                }
            }

//            if ( tGap(0,0) <= 0 && mResidualDofTypeRequested )
//            {
//                // compute the jacobian on slave side
//                // slave slave
//                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
//                        { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexSlave, 1 ) } )
//                        += tTempParam * trans( tFISlave -> N()) * tNormQ * tFISlave -> N() * tWStar;
//
//                // slave master
//                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
//                        { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tDofIndexMaster, 1 ) } )
//                        += (-1) * tTempParam * trans( tFISlave -> N()) * tNormQ * tFIMaster -> N() * tWStar;
//
//                // compute the jacobian on master side
//                //master master
//                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
//                        { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexMaster, 1 ) } )
//                        += tTempParam* trans( tFIMaster -> N()) * tNormQ * tFIMaster -> N() * tWStar;
//
//                //master slave
//                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
//                        { mSet->get_jac_dof_assembly_map()( tDofIndexMaster  )( tDofIndexSlave, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tDofIndexSlave, 1 ) } )
//                        += (-1) * tTempParam * trans( tFIMaster -> N()) * tNormQ * tFISlave -> N() * tWStar;
//            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_drdpdv( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Contact_Penalty::compute_drdpdv - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */


