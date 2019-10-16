
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

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
//			print(mMasterProp( 0 )->val(),"Neumann");

            // compute the residual r_U
//            aResidual( 0 ) = - trans( tShapeFunc2D ) * trans( mMasterProp( 0 )->val() );
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * trans( mMasterProp( 0 )->val() );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set jacobian size
            this->set_jacobian( aJacobians );

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

            // compute the jacobian for direct IWG dof dependencies
            // None

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency in the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
//                    aJacobians( 0 )( iDOF ).matrix_data()
//                    += - trans( tShapeFunc2D ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( 0 ) );
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                     moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
