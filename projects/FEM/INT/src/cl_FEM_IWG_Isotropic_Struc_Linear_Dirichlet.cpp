
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

        IWG_Isotropic_Struc_Linear_Dirichlet::IWG_Isotropic_Struc_Linear_Dirichlet()
        {

            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Dirichlet" ] = IWG_Property_Type::DIRICHLET;
            mPropertyMap[ "Select" ]    = IWG_Property_Type::SELECT;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "ElastLinIso" ] = IWG_Constitutive_Type::ELAST_LIN_ISO;

            // set size for the stabilization parameter pointer cell
            mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

            // populate the stabilization map
            mStabilizationMap[ "DirichletNitsche" ] = IWG_Stabilization_Type::DIRICHLET_NITSCHE;

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( mResidualDofType( 0 ) );

            // get SP, CM and property indices
            uint tSelectIndex      = static_cast< uint >( IWG_Property_Type::SELECT );
            uint tDirichletIndex   = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            uint tNitscheIndex     = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix if needed
            if ( mMasterProp( tSelectIndex ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();
                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( tSelectIndex )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

            // compute the residual
            mSet->get_residual()( 0 )( { tMasterResStartIndex, tMasterResStopIndex }, { 0, 0 } )
            += ( - trans( tFI->N() ) * tM * mMasterCM( tElastLinIsoIndex )->traction( mNormal )
                 + mMasterCM( tElastLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tM * tJump
                 + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * tJump ) * aWStar;

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get master index for residual dof type (here displacement), indices for assembly
            uint tMasterDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tMasterResStartIndex = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 0 );
            uint tMasterResStopIndex  = mSet->get_res_dof_assembly_map()( tMasterDofIndex )( 0, 1 );

            // get field interpolator for a given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get SP, CM and property indices
            uint tSelectIndex      = static_cast< uint >( IWG_Property_Type::SELECT );
            uint tDirichletIndex   = static_cast< uint >( IWG_Property_Type::DIRICHLET );
            uint tElastLinIsoIndex = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
            uint tNitscheIndex     = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );

            // selection matrix
            Matrix< DDRMat > tM;

            // set a default selection matrix
            if ( mMasterProp( tSelectIndex ) == nullptr )
            {
                // get spatial dimension
                uint tSpaceDim = tFI->get_dof_type().size();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = mMasterProp( tSelectIndex )->val();
            }

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( tDirichletIndex )->val();

//            // compute the jacobian for direct dof dependencies
//            if ( mResidualDofTypeRequested )
//            {
//                // compute the jacobian for direct dof dependencies
//                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
//                                      { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
//                += (   mMasterCM( tElastLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tM * tFI->N()
//                     + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * tFI->N() ) * aWStar;
//            }

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Master_Slave::MASTER );
                uint tMasterDepStartIndex = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 0 );
                uint tMasterDepStopIndex  = mSet->get_jac_dof_assembly_map()( tMasterDofIndex )( tDofDepIndex, 1 );

                // if dof type is residual dof type
                if( tDofType( 0 ) == mResidualDofType( 0 ) )
                {
                    // compute the jacobian for direct dof dependencies
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += (   mMasterCM( tElastLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tM * tFI->N()
                         + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * tFI->N() ) * aWStar;
                }

                // if dependency on the dof type
                if ( mMasterProp( tDirichletIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    -= ( mMasterCM( tElastLinIsoIndex )->testTraction( mNormal, mResidualDofType ) * tM * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType )
                       + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType )) * aWStar;
                }

                // if dependency on the dof type
                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                   -= aWStar * ( trans( tFI->N() ) *  tM * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) );
                }

                // if dependency on the dof type
                if ( mStabilizationParam( tNitscheIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { tMasterResStartIndex, tMasterResStopIndex },
                                          { tMasterDepStartIndex, tMasterDepStopIndex } )
                    += trans( tFI->N() ) * tM * tJump * mStabilizationParam( tNitscheIndex )->dSPdMasterDOF( tDofType ) * aWStar;
                }
            }

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
