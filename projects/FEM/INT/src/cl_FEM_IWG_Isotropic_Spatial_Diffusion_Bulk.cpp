
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            //FIXME heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // check master constitutive models
            this->check_constitutive_models();

            // compute flux
            Matrix< DDRMat > tStress;
            mMasterCM( 0 )->eval_stress( tStress );

            // set residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * tStress
                           - trans( mMasterFI( 0 )->N() ) * tQ;

//            // compute conductivity matrix
//            Matrix< DDRMat > I;
//            eye( mSpaceDim, mSpaceDim, I );
//            Matrix< DDRMat > K = mMasterProp( 0 )->val()( 0 ) * I;
//
//            // set residual size
//            this->set_residual( aResidual );
//
//            // compute the residual r_T
//            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->gradx( 1 )
//                           - trans( mMasterFI( 0 )->N() ) * tQ;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // check master constitutive models
            this->check_constitutive_models();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // compute the jacobian for residual dof type and direct IWG dof dependencies
                Matrix< DDRMat > tdStressdDOF;
                mMasterCM( 0 )->eval_dStressdDOF( mMasterGlobalDofTypes( iDOF ), tdStressdDOF );
                aJacobians( 0 )( iDOF ) = trans( mMasterFI( 0 )->Bx() ) * tdStressdDOF;
            }

//            // compute conductivity matrix
//            // FIXME in a constitutive model
//            Matrix< DDRMat > I;
//            eye( mSpaceDim, mSpaceDim, I );
//            Matrix< DDRMat > K = mMasterProp( 0 )->val()( 0 ) * I;
//
//            // set the jacobian size
//            this->set_jacobian( aJacobians );
//
//            // compute the jacobian for residual dof type and direct IWG dof dependencies
//            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->Bx();
//
//            // compute the jacobian for indirect IWG dof dependencies through properties
//            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
//            {
//                // if dependency in the dof type
//                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
//                {
//                    // compute property derivative
//                    // FIXME in a constitutive model
//                    Matrix< DDRMat > dK = mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( 0 ) );
//
//                    // add contribution to jacobian
//                    aJacobians( 0 )( iDOF ).matrix_data() += trans( mMasterFI( 0 )->Bx() ) * mMasterFI( 0 )->gradx( 1 ) * dK;
//                }
//            }
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            //FIXME heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // compute conductivity matrix
            // FIXME in a constitutive model
            Matrix< DDRMat > I;
            eye( mSpaceDim, mSpaceDim, I );
            Matrix< DDRMat > K = mMasterProp( 0 )->val()( 0 ) * I;

            // set the jacobian size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->gradx( 1 )
                           + trans( mMasterFI( 0 )->N() ) * tQ;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian j_T_T
            aJacobians( 0 )( 0 ) = trans( mMasterFI( 0 )->Bx() ) * K * mMasterFI( 0 )->Bx();

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // if dependency in the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute property derivative
                    // FIXME in a constitutive model
                    Matrix< DDRMat > dK = mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( 0 ) );

                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data() += trans( mMasterFI( 0 )->Bx() ) * mMasterFI( 0 )->gradx( 1 ) * dK;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
