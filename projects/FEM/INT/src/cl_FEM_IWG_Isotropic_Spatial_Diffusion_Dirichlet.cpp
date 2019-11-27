
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

        IWG_Isotropic_Spatial_Diffusion_Dirichlet::IWG_Isotropic_Spatial_Diffusion_Dirichlet()
        {
            // FIXME set a penalty
            mGamma = 1.0;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_residual( real tWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();
#endif

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( 0 )->val();

            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                          += ( - trans( tFI->N() ) * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mGamma * trans( tFI->N() ) * tJump )                * tWStar;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian( real tWStar )
        {
#ifdef DEBUG
            // check field interpolators, properties, constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();
#endif

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute jump
            Matrix< DDRMat > tJump = tFI->val() - mMasterProp( 0 )->val();

            if (mResidualDofTypeRequested)
            {
            // compute the jacobian for direct dof dependencies
            mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                  { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                                      += ( mMasterCM( 0 )->testTraction( mNormal ) * tFI->N()
                                        + mGamma * trans( tFI->N() ) * tFI->N() )                   * tWStar;
            }

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mRequestedMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mRequestedMasterGlobalDofTypes( iDOF );

                uint tIndexDep = mSet->get_dof_index_for_type( mRequestedMasterGlobalDofTypes( iDOF )( 0 ), mtk::Master_Slave::MASTER );

                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += ( -1.0 * mMasterCM( 0 )->testTraction( mNormal ) * mMasterProp( 0 )->dPropdDOF( tDofType )
                                 - mGamma * trans( tFI->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType ) )                       * tWStar;
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += ( - trans( tFI->N() ) * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal )
                                 + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) * tJump( 0 ) )           * tWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Dirichlet::compute_jacobian_and_residual - Not implemeted." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
