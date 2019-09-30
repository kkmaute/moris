
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

            // evaluate temperature jump
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

            // evaluate master conductivity matrix
            Matrix< DDRMat > tMasterK;
            mMasterCM( 0 )->eval_const( tMasterK );

            // evaluate slave conductivity matrix
            Matrix< DDRMat > tSlaveK;
            mSlaveCM( 0 )->eval_const( tSlaveK );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * mMasterFI( 0 )->N()
                                 + mGammaInterface * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 0 )( 1 ) = - trans( mMasterWeight * tMasterK * mMasterFI( 0 )->dnNdxn( 1 ) ) * mNormal * mSlaveFI( 0 )->N()
                                   - mGammaInterface * trans( mMasterFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            aJacobians( 1 )( 0 ) = trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * mMasterFI( 0 )->N()
                                   - mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 1 )( 1 ) = - trans( mSlaveWeight * tSlaveK * mSlaveFI( 0 )->dnNdxn( 1 ) ) * mNormal * mSlaveFI( 0 )->N()
                                   + mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            uint tMasterNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // evaluate master flux
                    Matrix< DDRMat > tMasterdStress;
                    mMasterCM( 0 )->eval_dStressdDOF( mMasterGlobalDofTypes( iDOF ), tMasterdStress );

                    // evaluate master contitutive matrix
                    Matrix< DDRMat > tMasterdConst;
                    mMasterCM( 0 )->eval_dConstdDOF( mMasterGlobalDofTypes( iDOF ), tMasterdConst );

                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterWeight * trans( mNormal ) * tMasterdStress
                           + mMasterWeight * trans( tMasterdConst ) * trans( mNormal ) * mMasterFI( 0 )->dnNdxn( 1 ) * tJump( 0 );
                    aJacobians( 1 )( iDOF ).matrix_data()
                        += trans( mSlaveFI( 0 )->N() ) * mMasterWeight * trans( mNormal ) * tMasterdStress;
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // if dependency on the dof type
                if ( mSlaveCM( 0 )->check_dof_dependency( mSlaveGlobalDofTypes( iDOF ) ) )
                {
                    // evaluate slave flux
                    Matrix< DDRMat > tSlavedStress;
                    mSlaveCM( 0 )->eval_dStressdDOF( mSlaveGlobalDofTypes( iDOF ), tSlavedStress );

                    // evaluate slave contitutive matrix
                    Matrix< DDRMat > tSlavedConst;
                    mSlaveCM( 0 )->eval_dConstdDOF( mSlaveGlobalDofTypes( iDOF ), tSlavedConst );

                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mSlaveWeight * trans( mNormal ) * tSlavedStress;
                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                        += trans( mSlaveFI( 0 )->N() ) * mSlaveWeight * trans( mNormal ) * tSlavedStress
                           + mSlaveWeight * trans( tSlavedConst )* trans( mNormal ) * mSlaveFI( 0 )->dnNdxn( 1 ) * tJump( 0 );
                }
            }

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
