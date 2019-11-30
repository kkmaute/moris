
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"

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
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();

            // compute jump
            Matrix< DDRMat > tJumpMat;
            this->build_jump( tJumpMat );
//            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - tJumpMat;

            Matrix<DDRMat> tConstDofs;
            this->get_I( tConstDofs );

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tConstDofs * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tConstDofs * tJump
                             + mGamma * trans( mMasterFI( 0 )->N() ) * tConstDofs * tJump;
        }
        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::build_jump( Matrix< DDRMat > & aJumpMat )
        {
            uint tDim = mMasterFI( 0 )->get_coeff().n_cols();

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
                        else if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UY )
                        {
                            aJumpMat(1) = mMasterProp( Ik )->val( )(0);
                        }
                        else if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UZ )
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
            uint tDim = mMasterFI( 0 )->get_coeff().n_cols();

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
                        else if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UY )
                        {
                            aI(1,1) = 1;
                        }
                        else if( mMasterProp( Ik )->get_dof_type( ) == MSI::Dof_Type::UZ )
                        {
                            aI(2,2) = 1;
                        }
                    }
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
//            this->check_properties();
//            this->check_constitutive_models();

            // compute jump
            Matrix< DDRMat > tJumpMat;
            this->build_jump( tJumpMat );

//            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - tJumpMat;

            Matrix<DDRMat> tConstDofs;
            this->get_I( tConstDofs );
            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( 0 )->testTraction( mNormal ) *tConstDofs* mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) *tConstDofs* mMasterFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through properties
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();

            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( 0 )->testTraction( mNormal ) *tConstDofs* mMasterProp( 0 )->dPropdDOF( tDofType )
                       - mGamma * trans( mMasterFI( 0 )->N() ) *tConstDofs* mMasterProp( 0 )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) *tConstDofs*  mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal ) ;
                      // + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump;
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
