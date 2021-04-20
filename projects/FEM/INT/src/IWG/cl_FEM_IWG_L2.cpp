#include "cl_FEM_IWG_L2.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "op_times.hpp" //LINALG/src
#include "fn_norm.hpp"  //LINALG/src
#include "fn_trans.hpp" //LINALG/src
#include "fn_dot.hpp"   //LINALG/src

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        IWG_L2::IWG_L2( const real aAlpha )
            {
            // set alpha
            this->set_alpha( aAlpha );
            }

        //------------------------------------------------------------------------------

        void IWG_L2::set_alpha( const real aAlpha )
        {
            mAlpha = aAlpha;

            if(  std::abs(aAlpha) < 0.0 )
            {
                mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_without_alpha;

                mComputeJacFunction
                = & IWG_L2::compute_jacobian_without_alpha;

                mComputeResFunction
                = & IWG_L2::compute_residual_without_alpha;
            }
            else
            {
                mComputeFunction
                = & IWG_L2::compute_jacobian_and_residual_with_alpha;

                mComputeJacFunction
                = & IWG_L2::compute_jacobian_with_alpha;

                mComputeResFunction
                = & IWG_L2::compute_residual_with_alpha;
            }
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual( real aWStar )
        {
            // call the right mComputeFunction for residual and jacobian evaluations
            ( this->*mComputeFunction )( aWStar );
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_without_alpha( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_jacobian_and_residual_without_alpha - will not work because of weights.");

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // if residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar * ( trans( tFI->N() ) * tFI->N() );

                    // compute residual
                    mSet->get_residual()(0)(
                            { tResStartIndex, tResStopIndex } ) +=
                                    mSet->get_jacobian()(
                                            { tResStartIndex, tResStopIndex },
                                            { tDepStartIndex, tDepStopIndex } ) * ( tFI->get_coeff() - mNodalWeakBCs );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_and_residual_with_alpha( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_L2::compute_jacobian_and_residual_with_alpha - will not work because of weights.");

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // if residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar *
                            ( trans( tFI->N() ) * tFI->N() + mAlpha * ( trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 ) ) );

                    // compute residual
                    mSet->get_residual()(0)(
                            { tResStartIndex, tResStopIndex } ) +=
                                    mSet->get_jacobian()(
                                            { tResStartIndex, tResStopIndex },
                                            { tDepStartIndex, tDepStopIndex } ) * ( tFI->get_coeff() - mNodalWeakBCs );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian( real aWStar )
        {
            // call the right mComputeJacFunction for jacobian evaluation
            ( this->*mComputeJacFunction )( aWStar );
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_without_alpha( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // if residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar * ( trans( tFI->N() ) * tFI->N() );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_jacobian_with_alpha( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // get the number of master dof type dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            // loop over master dof type dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofDepIndex, 1 );

                // if residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    // compute Jacobian
                    mSet->get_jacobian()(
                            { tResStartIndex, tResStopIndex },
                            { tDepStartIndex, tDepStopIndex } ) += aWStar *
                            ( trans( tFI->N() ) * tFI->N() + mAlpha * ( trans( tFI->dnNdxn( 1 ) ) * tFI->dnNdxn( 1 ) ) );
                }
            }
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_residual( real aWStar )
        {
            // call the right mComputeResFunction for residual evaluation
            ( this->*mComputeResFunction )( aWStar );
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_residual_without_alpha( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute residual
            mSet->get_residual()(0)(
                    { tResStartIndex, tResStopIndex } ) += aWStar *
                    trans( tFI->N() ) * ( tFI->val() - tFI->N() * mNodalWeakBCs ) ;
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_residual_with_alpha( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators
            this->check_field_interpolators();
#endif

            // get index for residual dof type, indices for assembly
            uint tDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartIndex = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tResStopIndex  = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // get the field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) ( 0 ));

            // compute residual
            mSet->get_residual()(0)(
                    { tResStartIndex, tResStopIndex } ) += aWStar *
                    ( trans( tFI->N() ) * ( tFI->val() - tFI->N() * mNodalWeakBCs ) +
                            mAlpha * trans( tFI->dnNdxn( 1 ) ) * ( tFI->gradx( 1 ) - tFI->dnNdxn( 1 ) * mNodalWeakBCs ) );
        }

        //------------------------------------------------------------------------------

        void IWG_L2::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_L2::compute_dRdp - not implemented." );
        }

        //------------------------------------------------------------------------------
        //
        //        real
        //        IWG_L2::interpolate_scalar_at_point(
        //                                    const Matrix< DDRMat > & aNodalWeakBC,
        //                                    const uint             & aPointIndex )
        //        {
        //            // get shape function
        //            mN->compute( aPointIndex );
        //
        //            // return interpolation
        //            return dot( mN->matrix() , aNodalWeakBC );
        //        }

        //------------------------------------------------------------------------------
        //
        //        real
        //        IWG_L2::compute_integration_error(
        //                const Matrix< DDRMat >                    & aNodalDOF,
        //                real (*aFunction)( const Matrix< DDRMat > & aPoint ) ,
        //                const uint                                & aPointIndex )
        //        {
        //            mN->compute( aPointIndex );
        //
        //            Matrix< DDRMat > tCoords = mN->matrix_data() * mInterpolator->get_node_coords();
        //
        //            // get shape function
        //            Matrix< DDRMat > tPhiHat = mN->matrix_data() * aNodalDOF.matrix_data();
        //
        //            return std::pow( tPhiHat( 0 ) - aFunction( tCoords ), 2 );
        //        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
