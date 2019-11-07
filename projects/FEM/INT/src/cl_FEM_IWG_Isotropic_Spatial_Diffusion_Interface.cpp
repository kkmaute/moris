
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
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

//            // check master and slave properties
//            this->check_properties( mtk::Master_Slave::MASTER );
//            this->check_properties( mtk::Master_Slave::SLAVE );
//
//            // check master and slave constitutive models
//            this->check_constitutive_models( mtk::Master_Slave::MASTER );
//            this->check_constitutive_models( mtk::Master_Slave::SLAVE );

            // set residual cell size
            this->set_residual_double( aResidual );

            Matrix< DDRMat > tTraction1 = mMasterCM( 0 )->traction( mNormal );

            Matrix< DDRMat > tTraction2 = mSlaveCM( 0 )->traction( mNormal );

            // evaluate average traction
            Matrix< DDRMat > tTraction = mMasterWeight * tTraction1 + mSlaveWeight * tTraction2;

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

//            print(mMasterCM( 0 )->traction( mNormal ),"mMasterCM( 0 )->traction( mNormal )");

            // compute master residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tTraction
                             + mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mGammaInterface * trans( mMasterFI( 0 )->N() ) * tJump;

            // compute slave residual
            aResidual( 1 ) =   trans( mSlaveFI( 0 )->N() ) * tTraction
                             + mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * tJump
                             - mGammaInterface * trans( mSlaveFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Interface::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

//            // check master and slave properties
//            this->check_properties( mtk::Master_Slave::MASTER );
//            this->check_properties( mtk::Master_Slave::SLAVE );
//
//            // check master and slave constitutive models
//            this->check_constitutive_models( mtk::Master_Slave::MASTER );
//            this->check_constitutive_models( mtk::Master_Slave::SLAVE );

            // set the jacobian cell size
            this->set_jacobian_double( aJacobians );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mMasterGlobalDofTypes.size();

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) =   mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                   + mGammaInterface * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 0 )( tMasterNumDofDependencies ) = - mMasterWeight * mMasterCM( 0 )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
                                                           - mGammaInterface * trans( mMasterFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            aJacobians( 1 )( 0 ) =   mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
                                   - mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 1 )( tMasterNumDofDependencies ) = - mSlaveWeight * mSlaveCM( 0 )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
                                                           + mGammaInterface * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            // compute the jacobian for indirect dof dependencies through master constitutive models
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // get the dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterWeight * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal )
                       + mMasterWeight * mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mMasterWeight * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal );
                }
            }

            // compute the jacobian for indirect dof dependencies through slave constitutive models
            uint tSlaveNumDofDependencies = mSlaveGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tSlaveNumDofDependencies; iDOF++ )
            {
                // get dof type
                Cell< MSI::Dof_Type > tDofType = mSlaveGlobalDofTypes( iDOF );

                // if dependency on the dof type
                if ( mSlaveCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mSlaveWeight * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mSlaveWeight * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal )
                       + mSlaveWeight * mSlaveCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump( 0 );
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
