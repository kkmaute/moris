
#include "cl_FEM_IWG_Diffusion_Dirichlet_Nitsche.hpp"
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

        IWG_Diffusion_Dirichlet_Nitsche::IWG_Diffusion_Dirichlet_Nitsche( sint aBeta )
        {
            // set sign for symmetric/unsymmetric Nitsche
            mBeta = aBeta;

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "Diffusion" ] = static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO );

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Dirichlet_Nitsche::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolatopr for a given dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the selection matrix property
            std::shared_ptr< Property > tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get imposed temperature property
            std::shared_ptr< Property > tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the elasticity CM
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - tPropDirichlet->val();

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tMasterResStartIndex, tMasterResStopIndex },
                    { 0, 0 } ) += aWStar * (
                            - trans( tFI->N() ) * tM * tCMDiffusion->traction( mNormal )
                            + mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType ) * tM * tJump
                            + tSPNitsche->val()( 0 ) * trans( tFI->N() ) * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type, indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator * tFI =
                    mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get the selection matrix property
            std::shared_ptr< Property > tPropSelect =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // get imposed temperature property
            std::shared_ptr< Property > tPropDirichlet =
                    mMasterProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get the elasticity CM
            std::shared_ptr< Constitutive_Model > tCMDiffusion =
                    mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) );

            // get the elasticity CM
            std::shared_ptr< Stabilization_Parameter > tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - tPropDirichlet->val();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > & tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    mBeta * tCMDiffusion->testTraction( mNormal, mResidualDofType ) * tM * tFI->N()
                                    + tSPNitsche->val()( 0 ) * trans( tFI->N() ) * tM * tFI->N() ) ;
                }

                // if dependency on the dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - mBeta  * tCMDiffusion->testTraction( mNormal, mResidualDofType ) * tM * tPropDirichlet->dPropdDOF( tDofType )
                                    - tSPNitsche->val()( 0 ) * trans( tFI->N() ) * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if dependency on the dof type
                if ( tCMDiffusion->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    - trans( tFI->N() ) * tM * tCMDiffusion->dTractiondDOF( tDofType, mNormal )
                                    + mBeta * tCMDiffusion->dTestTractiondDOF( tDofType, mNormal, mResidualDofType ) * tM( 0 ) * tJump( 0 ) );
                }

                // if dependency on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()(
                            { tMasterResStartIndex, tMasterResStopIndex },
                            { tMasterDepStartIndex, tMasterDepStopIndex } ) += aWStar * (
                                    trans( tFI->N() ) * tM * tJump( 0 ) * tSPNitsche->dSPdMasterDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Dirichlet_Nitsche::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void IWG_Diffusion_Dirichlet_Nitsche::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Dirichlet_Nitsche::compute_dRdp - Not implemented.");
        }

        //------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
