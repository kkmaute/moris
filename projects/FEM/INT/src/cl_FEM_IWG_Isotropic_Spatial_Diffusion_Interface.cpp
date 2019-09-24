
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Interface.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IWG_Isotropic_Spatial_Diffusion_Interface::IWG_Isotropic_Spatial_Diffusion_Interface()
        {
            // FIXME set penalty parameter
            mGammaInterface = 1.0;

            // FIXME set weight master parameter
            mMasterWeight = 0.5;

            // FIXME set weight slave parameter
            mSlaveWeight = 0.5;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );

            // check master and slave properties
            this->check_properties( mtk::Master_Slave::MASTER );
            this->check_properties( mtk::Master_Slave::SLAVE );

            // check master and slave constitutive models
            this->check_constitutive_models( mtk::Master_Slave::MASTER );
            this->check_constitutive_models( mtk::Master_Slave::SLAVE );

            // set residual cell size
            this->set_residual_double( aResidual );

            // evaluate master flux
            Matrix< DDRMat > tMasterStress;
            mMasterCM( 0 )->eval_stress( tMasterStress );

            // evaluate slave flux
            Matrix< DDRMat > tSlaveStress;
            mSlaveCM( 0 )->eval_stress( tSlaveStress );

            // evaluate average flux
            Matrix< DDRMat > tStress = mMasterWeight * tMasterStress + mSlaveWeight * tSlaveStress;

            // evaluate master conductivity matrix
            Matrix< DDRMat > tMasterK;
            mMasterCM( 0 )->eval_const( tMasterK );

            // evaluate slave conductivity matrix
            Matrix< DDRMat > tSlaveK;
            mSlaveCM( 0 )->eval_const( tSlaveK );

            // evaluate tempertaure jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute master residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * dot( tStress, mNormal )
                             +  trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * tJump
                             + mGammaInterface * trans( mMasterFI( 0 )->N() ) * tJump;

            // compute slave residual
            aResidual( 1 ) =   trans( mSlaveFI( 0 )->N() ) * dot( tStress, mNormal )
                             + trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * tJump
                             - mGammaInterface * trans( mSlaveFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master and slave field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_field_interpolators( mtk::Master_Slave::SLAVE );

            // check master and slave properties
            this->check_properties( mtk::Master_Slave::MASTER );
            this->check_properties( mtk::Master_Slave::SLAVE );

            // check master and slave constitutive models
            this->check_constitutive_models( mtk::Master_Slave::MASTER );
            this->check_constitutive_models( mtk::Master_Slave::SLAVE );

            // set the jacobian cell size
            this->set_jacobian_double( aJacobians );

            // evaluate master flux
            Matrix< DDRMat > tMasterdStress;
            mMasterCM( 0 )->eval_dStressdDOF( mResidualDofType, tMasterdStress );

            // evaluate slave flux
            Matrix< DDRMat > tSlavedStress;
            mSlaveCM( 0 )->eval_dStressdDOF( mResidualDofType, tSlavedStress );

            // evaluate master conductivity matrix
            Matrix< DDRMat > tMasterK;
            mMasterCM( 0 )->eval_const( tMasterK );

            // evaluate slave conductivity matrix
            Matrix< DDRMat > tSlaveK;
            mSlaveCM( 0 )->eval_const( tSlaveK );

            // evaluate master contitutive matrix
            Matrix< DDRMat > tMasterdConst;
            mMasterCM( 0 )->eval_dConstdDOF( mResidualDofType, tMasterdConst );

            // evaluate slave flux
            Matrix< DDRMat > tSlavedConst;
            mSlaveCM( 0 )->eval_dStressdDOF( mResidualDofType, tSlavedConst );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute Jacobian
            // FIXME add constitutive model dependencies
//            aJacobians( 0 )( 0 ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterWeight * trans( mNormal ) * tMasterdStress
//                                                  + mMasterWeight * tMasterdConst * trans( mNormal ) * mMasterFI( 0 )->dnNdxn( 1 ) * tJump
//                                                  + trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * mMasterFI( 0 )->N()
//                                                  + mGammaInterface * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();
            aJacobians( 0 )( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterWeight * trans( mNormal ) * tMasterdStress
                                                  + trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * mMasterFI( 0 )->N()
                                                  + mGammaInterface * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 0 )( 1 ) = - trans( mMasterFI( 0 )->N() ) * mSlaveWeight * trans( mNormal ) * tSlavedStress
                                                  - trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * mSlaveFI( 0 )->N()
                                                  - mGammaInterface * trans( mMasterFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            aJacobians( 1 )( 0 ) =   trans( mSlaveFI( 0 )->N() ) * mMasterWeight * trans( mNormal ) * tMasterdStress
                                                  + trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * mMasterFI( 0 )->N()
                                                  - mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mMasterFI( 0 )->N();

//            aJacobians( 1 )( 1 ).matrix_data() +=   trans( mSlaveFI( 0 )->N() ) * mSlaveWeight * trans( mNormal ) * tSlavedStress
//                                                  + mSlaveWeight * tSlavedConst * trans( mNormal ) * mSlaveFI( 0 )->dnNdxn( 1 ) * tJump
//                                                  - trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * mSlaveFI( 0 )->N()
//                                                  + mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();
            aJacobians( 1 )( 1 ) =   trans( mSlaveFI( 0 )->N() ) * mSlaveWeight * trans( mNormal ) * tSlavedStress
                                                  - trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * mSlaveFI( 0 )->N()
                                                  + mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();

        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                       moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
