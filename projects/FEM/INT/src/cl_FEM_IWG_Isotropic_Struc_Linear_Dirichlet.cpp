
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
            mGamma = 1.0;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

//            // compute flux
//            Matrix< DDRMat > tFlux;
//            mMasterCM( 0 )->eval_flux( tFlux );
//
//            // compute conductivity matrix
//            Matrix< DDRMat > tK;
//            mMasterCM( 0 )->eval_const( tK );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set residual size
            this->set_residual( aResidual );

//            Matrix< DDRMat > tNormal(2 ,3, 0.0 );
//            tNormal( 0, 0 ) = mNormal( 0,0 );
//            tNormal( 0, 2 ) = mNormal( 1,0 );
//            tNormal( 1, 1 ) = mNormal( 1,0 );
//            tNormal( 1, 2 ) = mNormal( 0,0 );
//
//            Matrix< DDRMat >tShapeFunc2D( 2, 8, 0.0 );
//            tShapeFunc2D( {0,0},{0,3} ) = mMasterFI( 0 )->N().matrix_data();
//            tShapeFunc2D( {1,1},{4,7} ) = mMasterFI( 0 )->N().matrix_data();
//
////			tShapeFunc2D( 0, 0 ) = mMasterFI( 0 )->N()( 0, 0 );
////            tShapeFunc2D( 1, 1 ) = mMasterFI( 0 )->N()( 0, 0 );
////            tShapeFunc2D( 0, 2 ) = mMasterFI( 0 )->N()( 0, 1 );
////            tShapeFunc2D( 1, 3 ) = mMasterFI( 0 )->N()( 0, 1 );
////            tShapeFunc2D( 0, 4 ) = mMasterFI( 0 )->N()( 0, 2 );
////            tShapeFunc2D( 1, 5 ) = mMasterFI( 0 )->N()( 0, 2 );
////            tShapeFunc2D( 0, 6 ) = mMasterFI( 0 )->N()( 0, 3 );
////            tShapeFunc2D( 1, 7 ) = mMasterFI( 0 )->N()( 0, 3 );
//
//            Matrix< DDRMat > tTestStrain;
//            mMasterCM( 0 )->eval_test_strain( tTestStrain );
//
//            // compute the residual
//            aResidual( 0 ) = - trans( tShapeFunc2D ) * tNormal * tFlux
//                             + trans( tTestStrain ) * tK * trans( tNormal ) * trans(tJump)
//                             + mGamma * trans( tShapeFunc2D ) * trans( tJump );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * trans( tJump )
                             + mGamma * trans( mMasterFI( 0 )->N() ) * trans( tJump );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();;

//            // compute conductivity matrix
//            Matrix< DDRMat > tK;
//            mMasterCM( 0 )->eval_const( tK );

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set the jacobian size
            this->set_jacobian( aJacobians );

//            Matrix< DDRMat > tNormal( 2, 3, 0.0 );
//            tNormal( 0, 0 ) = mNormal( 0,0 );
//            tNormal( 0, 2 ) = mNormal( 1,0 );
//            tNormal( 1, 1 ) = mNormal( 1,0 );
//            tNormal( 1, 2 ) = mNormal( 0,0 );
//
//            Matrix< DDRMat >tShapeFunc2D( 2, 8, 0.0 );
//            tShapeFunc2D( {0,0},{0,3} ) = mMasterFI( 0 )->N().matrix_data();
//            tShapeFunc2D( {1,1},{4,7} ) = mMasterFI( 0 )->N().matrix_data();
			
//			tShapeFunc2D( 0, 0 ) = mMasterFI( 0 )->N()( 0, 0 );
//            tShapeFunc2D( 1, 1 ) = mMasterFI( 0 )->N()( 0, 0 );
//            tShapeFunc2D( 0, 2 ) = mMasterFI( 0 )->N()( 0, 1 );
//            tShapeFunc2D( 1, 3 ) = mMasterFI( 0 )->N()( 0, 1 );
//            tShapeFunc2D( 0, 4 ) = mMasterFI( 0 )->N()( 0, 2 );
//            tShapeFunc2D( 1, 5 ) = mMasterFI( 0 )->N()( 0, 2 );
//            tShapeFunc2D( 0, 6 ) = mMasterFI( 0 )->N()( 0, 3 );
//            tShapeFunc2D( 1, 7 ) = mMasterFI( 0 )->N()( 0, 3 );
//
//            Matrix< DDRMat > tTestStrain;
//            mMasterCM( 0 )->eval_test_strain( tTestStrain );


            // compute the jacobian for direct dof dependencies
//            aJacobians( 0 )( 0 ) = trans( tTestStrain ) * tK * trans( tNormal ) * tShapeFunc2D
//                                   + mGamma * trans( tShapeFunc2D ) * tShapeFunc2D;
            print(mMasterFI( 0 )->N(),"mMasterFI( 0 )->N()");
            aJacobians( 0 )( 0 ) = mMasterCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                 + mGamma * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

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
//                    aJacobians( 0 )( iDOF ).matrix_data()
//                    += - trans( tTestStrain ) * tK * trans( tNormal ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) )
//                       - mGamma * trans( tShapeFunc2D ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += -1.0 * mMasterCM( 0 )->testTraction( mNormal ) * mMasterProp( 0 )->dPropdDOF( tDofType )
                       - mGamma * trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
//                    // evaluate stress derivative       //FIXME
//                    Matrix< DDRMat > tdFlux;
//                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdFlux );
//
//                    // evaluate constitutive matrix derivative
//                    Matrix< DDRMat > tdK;
//                    mMasterCM( 0 )->eval_dConstdDOF( mMasterGlobalDofTypes( iDOF ), tdK );

                    // add contribution to jacobian
//                    aJacobians( 0 )( iDOF ).matrix_data()
//                    += - trans( tShapeFunc2D ) * tNormal * tdFlux;
//                       //+ trans( tTestStrain ) * tNormal * tdK * tJump( 0 ) ;
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) *  mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal ) ;
                      // + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );
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
