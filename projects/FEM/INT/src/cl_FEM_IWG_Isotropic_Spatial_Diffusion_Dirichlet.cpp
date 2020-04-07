
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Dirichlet.hpp"
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
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( real tWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for a given dof type and indices for residual assembly
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartRow = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResEndRow   = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // get field interpolatopr for a given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get SP, CM, property indices
            uint tDirichletIndex  = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
            uint tNitscheIndex    = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

            // compute the residual
            mSet->get_residual()( 0 )( { tResStartRow, tResEndRow }, { 0, 0 } )
            += ( - trans( tFI->N() ) * mMasterCM( tDiffLinIsoIndex )->traction( mNormal )
                 + mMasterCM( tDiffLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tJump
                 + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tJump ) * tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get index for a given dof type and indices for residual assembly
            uint tResDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tResStartRow = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 0 );
            uint tResEndRow   = mSet->get_res_dof_assembly_map()( tResDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get SP, CM, property indices
            uint tDirichletIndex  = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tDiffLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );
            uint tNitscheIndex    = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

            // compute the jacobian for direct dof dependencies
            if ( mResidualDofTypeRequested )
            {
                mSet->get_jacobian()( { tResStartRow, tResEndRow },
                                      { mSet->get_jac_dof_assembly_map()( tResDofIndex )( tResDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tResDofIndex )( tResDofIndex, 1 ) } )
                += (   mMasterCM( tDiffLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tFI->N()
                     + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tFI->N() ) * tWStar;
            }

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex   = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tDepStartIndex = mSet->get_jac_dof_assembly_map()( tResDofIndex )( tDofDepIndex, 0 );
                uint tDepStopIndex  = mSet->get_jac_dof_assembly_map()( tResDofIndex )( tDofDepIndex, 1 );

                // if dependency on the dof type
                if ( mMasterProp( tDirichletIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tResStartRow,   tResEndRow }, { tDepStartIndex, tDepStopIndex } )
                    += ( -1.0 * mMasterCM( tDiffLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType )
                         - mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType ) ) * tWStar;
                }

                // if dependency on the dof type
                if ( mMasterCM( tDiffLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tResStartRow, tResEndRow }, { tDepStartIndex, tDepStopIndex } )
                    += ( - trans( tFI->N() ) * mMasterCM( tDiffLinIsoIndex )->dTractiondDOF( tDofType, mNormal )
                         + mMasterCM( tDiffLinIsoIndex )->dTestTractiondDOF( tDofType, mNormal, mResidualDofType ) * tJump( 0 ) ) * tWStar;
                }

                // if dependency on the dof type
                if ( mStabilizationParam( tNitscheIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tResStartRow, tResEndRow }, { tDepStartIndex, tDepStopIndex } )
                    += ( trans( tFI->N() ) * tJump( 0 ) * mStabilizationParam( tNitscheIndex )->dSPdMasterDOF( tDofType ) ) * tWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_dRdp - Not implemented.");
        }


//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
