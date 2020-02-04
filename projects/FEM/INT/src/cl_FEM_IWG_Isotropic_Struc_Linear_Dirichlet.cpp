
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
            mPropertyMap[ "Select" ] = IWG_Property_Type::SELECT;

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

            // get index for given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // get field interpolator for given dof type
            Field_Interpolator * tFI = mMasterFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

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

            // get start and end indices for residual assembly
            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
            += ( - trans( tFI->N() ) * tM * mMasterCM( tElastLinIsoIndex )->traction( mNormal )
                 + mMasterCM( tElastLinIsoIndex )->testTraction( mNormal ) * tM * tJump
                 + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * tJump ) * aWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef DEBUG
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get index for a given dof type
            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

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

            // compute the jacobian for direct dof dependencies
            if ( mResidualDofTypeRequested )
            {
                // compute the jacobian for direct dof dependencies
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                      { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                += (   mMasterCM( tElastLinIsoIndex )->testTraction( mNormal ) * tM * tFI->N()
                     + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * tFI->N() ) * aWStar;
            }

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                // get the index for this dof type
                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if dependency on the dof type
                if ( mMasterProp( tDirichletIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                    -= (mMasterCM( tElastLinIsoIndex )->testTraction( mNormal ) * tM * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType )
                       + mStabilizationParam( tNitscheIndex )->val()( 0 ) * trans( tFI->N() ) * tM * mMasterProp( tDirichletIndex )->dPropdDOF( tDofType )) * aWStar;
                }

                // if dependency on the dof type
                if ( mMasterCM( tElastLinIsoIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                   += ( - trans( tFI->N() ) *  tM * mMasterCM( tElastLinIsoIndex )->dTractiondDOF( tDofType, mNormal ) ) * aWStar;
                }

                // if dependency on the dof type
                if ( mStabilizationParam( tNitscheIndex )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
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
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_drdpdv( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_drdpdv - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
