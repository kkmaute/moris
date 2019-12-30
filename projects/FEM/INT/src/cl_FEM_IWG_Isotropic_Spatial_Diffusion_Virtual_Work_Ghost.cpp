
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

            // FIXME set penalty from outside
            // penalty parameter
            mGammaGhost = 1.0;

            // FIXME set mesh parameter from outside
            // mesh parameter
            mMeshParameter = 1.0;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

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

            // check, if order is supported
            MORIS_ERROR( mOrder < 4, " IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_residual - order not supported. " );

            // get indices for SP, CM and properties
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            if( mOrder >= 1 )
            {
                // penalty parameter
                real tGhostPenalty1 = mGammaGhost * mMeshParameter;

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( 1 );

                // compute the jump in traction
                Matrix< DDRMat > tGradJump = mMasterCM( tDiffLinIsoIndex )->traction( mNormal ) - mSlaveCM( tDiffLinIsoIndex )->traction( mNormal );

                // compute the residual
                mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) }, { 0, 0 } )
                +=   tGhostPenalty1 * trans( tMasterFI->dnNdxn( 1 ) ) * trans( tNormalMatrix ) * tGradJump * aWStar;

                mSet->get_residual()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) }, { 0, 0 } )
                += - tGhostPenalty1 * trans( tSlaveFI->dnNdxn( 1 ) )  * trans( tNormalMatrix ) * tGradJump * aWStar;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );
#endif

            // get master index for residual dof type
            uint tDofIndexMaster = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get slave index for residual dof type
            uint tDofIndexSlave  = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::SLAVE );

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

            // check, if order is supported
            MORIS_ERROR( mOrder < 4, " IWG_Isotropic_Spatial_Diffusion_Ghost::compute_jacobian - Ghost stabilization for this order not supported yet. " );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // get indices for SP, CM and properties
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // order 1
            if ( mOrder >= 1 )
            {
                // get normal matrix
                Matrix< DDRMat > tNormalMatrix = this->get_normal_matrix( 1 );

                // penalty parameter
                real tGhostPenalty = mGammaGhost * mMeshParameter;

                // compute the jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );

                    // if dependency on the dof type
                    if ( mMasterCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
                        +=   tGhostPenalty * trans( tMasterFI->dnNdxn( 1 ) ) * trans( tNormalMatrix ) * mMasterCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 ) } )
                        += - tGhostPenalty * trans( tSlaveFI->dnNdxn( 1 ) ) * trans( tNormalMatrix ) * mMasterCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * aWStar;
                    }
                }

                // compute the jacobian for indirect dof dependencies through slave constitutive models
                uint tSlaveNumDofDependencies = mRequestedSlaveGlobalDofTypes.size();
                for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
                {
                    // get dof type
                    Cell< MSI::Dof_Type > tDofType = mRequestedSlaveGlobalDofTypes( iDOF );

                    // get index for the dof type
                    sint tIndexDep = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::SLAVE );

                    // if dependency on the dof type
                    if ( mSlaveCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexMaster )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
                        += - tGhostPenalty * trans( tMasterFI->dnNdxn( 1 ) ) * trans( tNormalMatrix ) * mSlaveCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * aWStar;

                        mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndexSlave )( 0, 1 ) },
                                              { mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndexSlave )( tIndexDep, 1 ) } )
                        +=   tGhostPenalty * trans( tSlaveFI->dnNdxn( 1 ) ) * trans( tNormalMatrix ) * mSlaveCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) * aWStar;
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
        void IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                                moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
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
