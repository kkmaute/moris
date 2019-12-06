
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
//------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Dirichlet::IWG_Isotropic_Struc_Linear_Dirichlet()
        {
            // FIXME set a penalty
            mGamma = 1000000.0;

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( real aWStar )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute jump
            Matrix< DDRMat > tJumpMat;
            this->build_jump( tJumpMat );
//            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();
            Matrix< DDRMat > tJump = tFI->val() - tJumpMat;

            Matrix<DDRMat> tConstDofs;
            this->get_I( tConstDofs );

            uint tStartRow = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 );
            uint tEndRow   = mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 );

            // compute the residual
            mSet->get_residual()( { tStartRow, tEndRow }, { 0, 0 } )
                    += ( - trans( tFI->N() ) * tConstDofs * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tConstDofs * tJump
                             + mGamma * trans( tFI->N() ) * tConstDofs * tJump ) * aWStar;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Dirichlet::build_jump( Matrix< DDRMat > & aJumpMat )
        {
            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tDim = tFI->get_coeff().n_cols();

            aJumpMat.set_size( tDim, 1, 0.0 );

            for( uint Ik=0; Ik<mMasterProp.size(); ++Ik )
            {
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = this->get_dof_type_list();

                for( uint Ii=0; Ii<tDofTypes.size(); ++Ii )
                {
                    for( uint Ij=0; Ij<tDofTypes.size(); ++Ij )
                    {
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UX )
                        {
                            aJumpMat(0) = mMasterProp( Ik )->val( )(0);
                        }
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UY )
                        {
                            aJumpMat(1) = mMasterProp( Ik )->val( )(0);
                        }
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UZ )
                        {
                            aJumpMat(2) = mMasterProp( Ik )->val( )(0);
                        }
                    }

                }
            }
        }
         //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::get_I( Matrix< DDRMat > & aI )
        {
            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );
            uint tDim = tFI->get_coeff().n_cols();

            aI.set_size( tDim, tDim, 0.0 );

            for( uint Ik=0; Ik<mMasterProp.size(); ++Ik )
            {
                moris::Cell< moris::Cell< MSI::Dof_Type > > tDofTypes = this->get_dof_type_list();

                for( uint Ii=0; Ii<tDofTypes.size(); ++Ii )
                {
                    for( uint Ij=0; Ij<tDofTypes.size(); ++Ij )
                    {
//                        if( mMasterProp( Ik )->get_dof_type( ) == tDofTypes( Ii )( Ij ) )
//                        {
//                            tI( Ij,Ij ) = 1;
//                        }
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UX )
                        {
                            aI(0,0) = 1;
                        }
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UY )
                        {
                            aI(1,1) = 1;
                        }
                        if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UZ )
                        {
                            aI(2,2) = 1;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( real aWStar )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();

            Field_Interpolator * tFI = mFieldInterpolatorManager->get_field_interpolators_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            // compute jump
            Matrix< DDRMat > tJumpMat;
            this->build_jump( tJumpMat );

//            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();
            Matrix< DDRMat > tJump = tFI->val() - tJumpMat;

            Matrix<DDRMat> tConstDofs;
            this->get_I( tConstDofs );

            uint tDofIndex = mSet->get_dof_index_for_type( mResidualDofType( 0 ), mtk::Master_Slave::MASTER );

            if (mResidualDofTypeRequested)
            {
                // compute the jacobian for direct dof dependencies
                mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                      { mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tDofIndex, 1 ) } )
                                    += ( mMasterCM( 0 )->testTraction( mNormal ) *tConstDofs* tFI->N()
                                        + mGamma * trans( tFI->N() ) *tConstDofs* tFI->N() )              * aWStar;
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
                            += ( -1.0 * mMasterCM( 0 )->testTraction( mNormal ) *tConstDofs* mMasterProp( 0 )->dPropdDOF( tDofType )
                                 - mGamma * trans( tFI->N() ) *tConstDofs* mMasterProp( 0 )->dPropdDOF( tDofType ) )              * aWStar;
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    mSet->get_jacobian()( { mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 0 ), mSet->get_res_dof_assembly_map()( tDofIndex )( 0, 1 ) },
                                          { mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 0 ), mSet->get_jac_dof_assembly_map()( tDofIndex )( tIndexDep, 1 ) } )
                            += ( - trans( tFI->N() ) *tConstDofs*  mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal )
                               + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal, tJump ) )              * aWStar;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
