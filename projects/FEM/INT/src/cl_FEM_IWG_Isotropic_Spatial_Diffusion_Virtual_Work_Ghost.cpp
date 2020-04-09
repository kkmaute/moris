
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost.hpp"
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
        IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost()
        {
            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize(  static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostVWOrder" ] = IWG_Stabilization_Type::GHOST_VW;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tDofIndexMaster      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tDofIndexSlave      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tMasterFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tSlaveFI  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // FIXME the order should be set differently
            switch ( tMasterFI->get_space_interpolation_order() )
            {
                case( mtk::Interpolation_Order::LINEAR ):
                {
                    mOrder = 1;
                    break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ):
                {
                    mOrder = 2;
                    break;
                }
                case( mtk::Interpolation_Order::CUBIC ):
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // get the diffusion constitutive model for master and slave
            std::shared_ptr< Constitutive_Model > tCMMasterDiff
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveDiff
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPGhost
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for master and slave
            real tMasterConductivity = tCMMasterDiff->constitutive()( 0, 0 );
            real tSlaveConductivity  = tCMSlaveDiff->constitutive()( 0, 0 );

            // loop over the interpolation order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                 // get normal matrix
                 Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );

                 // multiply common terms
                 Matrix< DDRMat > tPreMultiply
                 = tSPGhost->val()( 0 ) * trans( tNormalMatrix ) * tNormalMatrix
                 * ( tMasterConductivity * tMasterFI->gradx( iOrder ) - tSlaveConductivity * tSlaveFI->gradx( iOrder ) );

                 // compute master residual
                 mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                 += aWStar * ( trans( tMasterFI->dnNdxn( iOrder ) ) * tPreMultiply );

                 // compute slave residual
                 mSet->get_residual()( 0 )( { tSlaveResStartIndex, tSlaveResStopIndex }, { 0, 0 } )
                 -= aWStar * ( trans( tSlaveFI->dnNdxn( iOrder ) ) * tPreMultiply );
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type, indices for assembly
            uint tDofIndexMaster      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 );

            // get slave index for residual dof type, indices for assembly
            uint tDofIndexSlave      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );
            uint tSlaveResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 );
            uint tSlaveResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 );

            // get master field interpolator for the residual dof type
            Field_Interpolator * tMasterFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get slave field interpolator for the residual dof type
            Field_Interpolator * tSlaveFI  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // FIXME the order should be set differently
            switch ( tMasterFI->get_space_interpolation_order() )
            {
                case( mtk::Interpolation_Order::LINEAR ):
                {
                    mOrder = 1;
                    break;
                }
                case( mtk::Interpolation_Order::QUADRATIC ):
                {
                    mOrder = 2;
                    break;
                }
                case( mtk::Interpolation_Order::CUBIC ):
                {
                    mOrder = 3;
                    break;
                }
                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual - order not supported");
                    break;
                }
            }

            // get the diffusion constitutive model for master and slave
            std::shared_ptr< Constitutive_Model > tCMMasterDiff
            = mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );
            std::shared_ptr< Constitutive_Model > tCMSlaveDiff
            = mSlaveCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the stabilization parameter
            std::shared_ptr< Stabilization_Parameter > tSPGhost
            = mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // FIXME get the conductivity for master and slave
            real tMasterConductivity = tCMMasterDiff->constitutive()( 0, 0 );
            real tSlaveConductivity  = tCMSlaveDiff->constitutive()( 0, 0 );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // get normal matrix
                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );

                // set the order for the stabilization parameter
                tSPGhost->set_interpolation_order( iOrder );

                // compute the jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        += aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tMasterConductivity * tMasterFI->dnNdxn( iOrder ) );

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        -= aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tMasterConductivity * tMasterFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMMasterDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        -= aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tMasterFI->gradx( iOrder ) * tCMMasterDiff->dConstdDOF( tDofType ) );

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        += aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tMasterFI->gradx( iOrder ) * tCMMasterDiff->dConstdDOF( tDofType ) );
                    }
                }

                // compute the jacobian for indirect dof dependencies through slave constitutive models
                uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                    // get index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 );

                    // if dependency on the dof type
                    if ( tDofType( 0 ) == mResidualDofType( 0 ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        -= aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tSlaveConductivity * tSlaveFI->dnNdxn( iOrder ) );

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        += aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tSlaveConductivity * tSlaveFI->dnNdxn( iOrder ) );
                    }

                    // if diffusion CM depends on dof type
                    if( tCMSlaveDiff->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        -= aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tMasterFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tSlaveFI->gradx( iOrder ) * tCMSlaveDiff->dConstdDOF( tDofType ) );

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        += aWStar * (   tSPGhost->val()( 0 )
                                      * trans( tSlaveFI->dnNdxn( iOrder ) ) * trans( tNormalMatrix )
                                      * tNormalMatrix * tSlaveFI->gradx( iOrder ) * tCMSlaveDiff->dConstdDOF( tDofType ) );
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_dRdp - This function does nothing.");
        }

//------------------------------------------------------------------------------
        Matrix< DDRMat > IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::get_normal_matrix ( uint aOrderGhost )
        {
            // init the normal matrix
            Matrix< DDRMat > tNormalMatrix;

            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // switch on the ghost order
            switch( aOrderGhost )
            {
                case ( 1 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            tNormalMatrix = trans( mNormal );
                            break;
                        }
                        case ( 3 ):
                        {
                            tNormalMatrix = trans( mNormal );
                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case ( 2 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 2, 3, 0.0 );

                            // fill the normal matrix
                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );

                            tNormalMatrix( 0, 2 ) = mNormal( 1 );
                            tNormalMatrix( 1, 2 ) = mNormal( 0 );

                            break;
                        }
                        case ( 3 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 3, 6, 0.0 );

                            // fill the normal matrix
                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );
                            tNormalMatrix( 2, 2 ) = mNormal( 2 );

                            tNormalMatrix( 1, 3 ) = mNormal( 2 );
                            tNormalMatrix( 2, 3 ) = mNormal( 1 );

                            tNormalMatrix( 0, 4 ) = mNormal( 2 );
                            tNormalMatrix( 2, 4 ) = mNormal( 0 );

                            tNormalMatrix( 0, 5 ) = mNormal( 1 );
                            tNormalMatrix( 1, 5 ) = mNormal( 0 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                case ( 3 ):
                {
                    switch ( tSpaceDim )
                    {
                        case ( 2 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 3, 4, 0.0 );

                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );

                            tNormalMatrix( 0, 2 ) = mNormal( 1 );
                            tNormalMatrix( 1, 3 ) = mNormal( 0 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            tNormalMatrix( 2, 2 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 2, 3 ) = tSqrtOf2 * mNormal( 1 );

                            break;
                        }
                        case ( 3 ):
                        {
                            // set the normal matrix size
                            tNormalMatrix.set_size( 6, 10, 0.0 );

                            tNormalMatrix( 0, 0 ) = mNormal( 0 );
                            tNormalMatrix( 1, 1 ) = mNormal( 1 );
                            tNormalMatrix( 2, 2 ) = mNormal( 2 );

                            tNormalMatrix( 0, 3 ) = mNormal( 1 );
                            tNormalMatrix( 0, 4 ) = mNormal( 2 );

                            tNormalMatrix( 1, 5 ) = mNormal( 0 );
                            tNormalMatrix( 1, 6 ) = mNormal( 2 );

                            tNormalMatrix( 2, 7 ) = mNormal( 0 );
                            tNormalMatrix( 2, 8 ) = mNormal( 1 );

                            real tSqrtOf2 = std::sqrt( 2 );

                            tNormalMatrix( 3, 3 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 3, 5 ) = tSqrtOf2 * mNormal( 1 );
                            tNormalMatrix( 3, 9 ) = tSqrtOf2 * mNormal( 2 );

                            tNormalMatrix( 4, 6 ) = tSqrtOf2 * mNormal( 1 );
                            tNormalMatrix( 4, 8 ) = tSqrtOf2 * mNormal( 2 );
                            tNormalMatrix( 4, 9 ) = tSqrtOf2 * mNormal( 0 );

                            tNormalMatrix( 5, 4 ) = tSqrtOf2 * mNormal( 0 );
                            tNormalMatrix( 5, 7 ) = tSqrtOf2 * mNormal( 2 );
                            tNormalMatrix( 5, 9 ) = tSqrtOf2 * mNormal( 1 );

                            break;
                        }
                        default:
                        {
                            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::get_normal_matrix - Spatial dimensions can only be 2, 3." );
                            break;
                        }
                    }
                    break;
                }

                default:
                {
                    MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::get_normal_matrix - order not supported." );
                    break;
                }
            }
            return tNormalMatrix;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
