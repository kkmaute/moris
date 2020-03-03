
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
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostVWOrder1" ] = IWG_Stabilization_Type::GHOST_VW_1;
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

//            // get slave field interpolator for the residual dof type
//            Field_Interpolator * tSlaveFI  = mSlaveFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

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
            MORIS_ERROR( mOrder <= 1, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost:compute_residual - only first order supported. ");

            // get indices for SP, CM and properties
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // penalty parameter
                uint tGhostPenalty = mStabilizationParam( iOrder - 1 )->val()( 0 );

                // get flattened normal matrix
                Matrix<DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );

                // compute the jump in traction
                Matrix< DDRMat > tGradJump = mMasterCM( tDiffLinIsoIndex )->traction( mNormal )
                                           - mSlaveCM( tDiffLinIsoIndex )->traction( mNormal );

                // compute the master residual
                mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
                += tGhostPenalty * trans( mMasterCM( tDiffLinIsoIndex )->testStrain() ) * trans( tNormalMatrix ) * tGradJump * aWStar;

                // compute the slave residual
                mSet->get_residual()( 0 )( { tSlaveResStartIndex,  tSlaveResStopIndex },  { 0, 0 } )
                -= tGhostPenalty * trans( mSlaveCM( tDiffLinIsoIndex )->testStrain() )  * trans( tNormalMatrix ) * tGradJump * aWStar;
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
            MORIS_ERROR( mOrder <= 1, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost:compute_residual - only first order supported. ");

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // get indices for SP, CM and properties
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // get normal matrix
                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( 1 );

                // penalty parameter
                real tGhostPenalty = mStabilizationParam( iOrder - 1 )->val()( 0 );

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
                    if ( mMasterCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        +=   tGhostPenalty * trans( tMasterFI->dnNdxn( 1 ) ) * trans( tNormalMatrix )
                                           * mMasterCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal )
                                           * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        += - tGhostPenalty * trans( tSlaveFI->dnNdxn( 1 ) ) * trans( tNormalMatrix )
                                           * mMasterCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal )
                                           * aWStar;
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
                    if ( mSlaveCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                              { tDepStartIndex,       tDepStopIndex } )
                        += - tGhostPenalty * trans( tMasterFI->dnNdxn( 1 ) ) * trans( tNormalMatrix )
                                           * mSlaveCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal )
                                           * aWStar;

                        mSet->get_jacobian()( { tSlaveResStartIndex, tSlaveResStopIndex },
                                              { tDepStartIndex,      tDepStopIndex } )
                        +=   tGhostPenalty * trans( tSlaveFI->dnNdxn( 1 ) ) * trans( tNormalMatrix )
                                           * mSlaveCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal )
                                           * aWStar;
                    }
                }
            }

            // FIXME higher orders
//            // loop over the interpolation orders
//            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
//            {
//                // get normal matrix
//                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( iOrder );
//
//                // penalty parameter
//                real tGhostPenalty = mGammaGhost * std::pow( mMeshParameter, 2 * ( iOrder - 1 ) + 1 );
//
//                // compute Jacobian direct dependencies
//                aJacobians( 0 )( 0 ).matrix_data()
//                +=   tGhostPenalty * trans( mMasterFI( 0 )->dnNdxn( iOrder ) ) * trans( tNormalMatrix ) * mMasterCM( 0 )->constitutive() * tNormalMatrix * mMasterFI( 0 )->dnNdxn( iOrder );
//                aJacobians( 0 )( tMasterNumDofDependencies ).matrix_data()
//                += - tGhostPenalty * trans( mMasterFI( 0 )->dnNdxn( iOrder ) ) * trans( tNormalMatrix ) * mSlaveCM( 0 )->constitutive() * tNormalMatrix * mSlaveFI( 0 )->dnNdxn( iOrder );
//
//                aJacobians( 1 )( 0 ).matrix_data()
//                += - tGhostPenalty * trans( mSlaveFI( 0 )->dnNdxn( iOrder ) ) * trans( tNormalMatrix ) * mMasterCM( 0 )->constitutive() * tNormalMatrix * mMasterFI( 0 )->dnNdxn( iOrder );
//                aJacobians( 1 )( tMasterNumDofDependencies ).matrix_data()
//                +=   tGhostPenalty * trans( mSlaveFI( 0 )->dnNdxn( iOrder ) ) * trans( tNormalMatrix ) * mSlaveCM( 0 )->constitutive() * tNormalMatrix * mSlaveFI( 0 )->dnNdxn( iOrder );
//
//                // compute the jacobian for indirect dof dependencies through master constitutive models
//                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
//                {
//                    // get the dof type
//                    Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );
//
//                    // if dependency on the dof type
//                    if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
//                    {
//
//                        Matrix< DDRMat > tPreMultiply = tGhostPenalty * trans( tNormalMatrix ) * tNormalMatrix;
//
//                        // add contribution to jacobian
//                        aJacobians( 0 )( iDOF ).matrix_data()
//                        +=   trans( mMasterFI( 0 )->dnNdxn( iOrder ) ) * tPreMultiply * mMasterFI( 0 )->gradx( iOrder ) * mMasterCM( 0 )->dConstdDOF( tDofType );
//
//                        aJacobians( 1 )( iDOF ).matrix_data()
//                        += - trans( mSlaveFI( 0 )->dnNdxn( iOrder ) ) * tPreMultiply * mMasterFI( 0 )->gradx( iOrder ) * mMasterCM( 0 )->dConstdDOF( tDofType );
//                    }
//                }
//
//                // compute the jacobian for indirect dof dependencies through slave constitutive models
//                uint tSlaveNumDofDependencies = mSlaveGlobalDofTypes.size();
//                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
//                {
//                    // get dof type
//                    Cell< MSI::Dof_Type > tDofType = mSlaveGlobalDofTypes( iDOF );
//
//                    // if dependency on the dof type
//                    if ( mSlaveCM( 0 )->check_dof_dependency( tDofType ) )
//                    {
//
//                        Matrix< DDRMat > tPreMultiply = tGhostPenalty * trans( tNormalMatrix ) * tNormalMatrix;
//
//                        // add contribution to jacobian
//                        aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
//                        += - trans( mMasterFI( 0 )->dnNdxn( iOrder ) ) * tPreMultiply * mSlaveFI( 0 )->gradx( iOrder ) * mSlaveCM( 0 )->dConstdDOF( tDofType );
//
//                        aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
//                        +=   trans( mSlaveFI( 0 )->dnNdxn( iOrder ) ) * tPreMultiply * mSlaveFI( 0 )->gradx( iOrder ) * mSlaveCM( 0 )->dConstdDOF( tDofType );
//                    }
//                }
//            }
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
