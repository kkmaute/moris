
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost.hpp"
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

        IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost()
        {
            // set ghost flag
            mIsGhost = true;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );
            mSlaveCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "GhostVW" ] = static_cast< uint >( IWG_Stabilization_Type::GHOST_VW );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_residual( real aWStar )
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

            // get indices for SP, CM and properties
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // loop over the order
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // FIXME check to be removed
                MORIS_ERROR( iOrder <= 1, "Not implemented for order higher than one yet" );

                // set the order for the stabilization parameter
                tSP->set_interpolation_order( iOrder );

                // get flattened normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, 2 );

                // compute the jump in traction
                Matrix< DDRMat > tGradJump =
                        mMasterCM( tElastLinIsoIndex )->traction( mNormal ) -
                        mSlaveCM( tElastLinIsoIndex )->traction( mNormal );

                // compute the residual
                mSet->get_residual()( 0 )(
                        { tMasterResStartIndex, tMasterResStopIndex },
                        { 0, 0 } ) += aWStar * (
                                tSP->val()( 0 ) *
                                trans( mMasterCM( tElastLinIsoIndex )->testStrain() ) *
                                trans( tNormalMatrix ) * tGradJump );

                mSet->get_residual()( 0 )(
                        { tSlaveResStartIndex,  tSlaveResStopIndex },
                        { 0, 0 } ) -= aWStar * (
                                tSP->val()( 0 ) *
                                trans( mSlaveCM( tElastLinIsoIndex )->testStrain() ) *
                                trans( tNormalMatrix ) * tGradJump );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian( real aWStar )
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

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // get indices for SP, CM and properties
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

            // get the stabilization parameter
            const std::shared_ptr< Stabilization_Parameter > & tSP =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::GHOST_VW ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // order 1
            for ( uint iOrder = 1; iOrder <= mOrder; iOrder++ )
            {
                // FIXME check to be removed
                MORIS_ERROR( iOrder <= 1, "Not implemented for order higher than one yet" );

                // set the order for the stabilization parameter
                tSP->set_interpolation_order( iOrder );

                // get normal matrix
                Matrix< DDRMat > tNormalMatrix;
                this->get_flat_normal_matrix( tNormalMatrix, 2 );

                // compute the jacobian for indirect dof dependencies through master constitutive models
                for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
                {
                    // get the dof type
                    const Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                    // get the index for the dof type
                    sint tIndexDep      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                    uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 0 );
                    uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndexMaster )( tIndexDep, 1 );

                    // if dependency on the dof type
                    if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tDepStartIndex,       tDepStopIndex } ) += aWStar * (
                                        tSP->val()( 0 ) *
                                        trans( mMasterCM( tElastLinIsoIndex )->testStrain() ) * trans( tNormalMatrix ) *
                                        mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tDepStartIndex,      tDepStopIndex } ) -= aWStar * (
                                        tSP->val()( 0 ) *
                                        trans( mSlaveCM( tElastLinIsoIndex )->testStrain() ) * trans( tNormalMatrix ) *
                                        mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) );
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
                    if ( mSlaveCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to jacobian
                        mSet->get_jacobian()(
                                { tMasterResStartIndex, tMasterResStopIndex },
                                { tDepStartIndex,       tDepStopIndex } ) -= aWStar * (
                                        tSP->val()( 0 ) *
                                        trans( mMasterCM( tElastLinIsoIndex )->testStrain() ) * trans( tNormalMatrix ) *
                                        mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) );

                        mSet->get_jacobian()(
                                { tSlaveResStartIndex, tSlaveResStopIndex },
                                { tDepStartIndex,      tDepStopIndex } ) += aWStar * (
                                        tSP->val()( 0 ) *
                                        trans( mSlaveCM( tElastLinIsoIndex )->testStrain() ) * trans( tNormalMatrix ) *
                                        mSlaveCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Virtual_Work_Ghost::get_flat_normal_matrix(
                       Matrix< DDRMat > & aFlatNormal,
                       uint               aOrder )
        {
            // get spatial dimensions
            uint tSpaceDim = mNormal.numel();

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

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
                            aFlatNormal.set_size( 2, tCMElasticity->strain().numel(), 0.0 );

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
