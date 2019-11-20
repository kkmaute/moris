
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

            // set residual cell size
            this->set_residual_double( aResidual );

            // evaluate average traction
            Matrix< DDRMat > tTraction = mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->traction( mNormal )
                                       + mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->traction( mNormal );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute master residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tTraction
                             + mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->testTraction( mNormal ) * tJump
                             + mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tJump;

            // compute slave residual
            aResidual( 1 ) =   trans( mSlaveFI( 0 )->N() ) * tTraction
                             + mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->testTraction( mNormal ) * tJump
                             - mStabilizationParam( 0 )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master and slave field interpolators
            this->check_dof_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dof_field_interpolators( mtk::Master_Slave::SLAVE );
            this->check_dv_field_interpolators( mtk::Master_Slave::MASTER );
            this->check_dv_field_interpolators( mtk::Master_Slave::SLAVE );

            // set the jacobian cell size
            this->set_jacobian_double( aJacobians );

            // get number of master dof dependencies
            uint tMasterNumDofDependencies = mMasterGlobalDofTypes.size();

            // evaluate displacement jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 )
            =   mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
              + mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 0 )( tMasterNumDofDependencies )
            = - mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
              - mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * mSlaveFI( 0 )->N();

            aJacobians( 1 )( 0 )
            =   mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->testTraction( mNormal ) * mMasterFI( 0 )->N()
              - mStabilizationParam( 0 )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * mMasterFI( 0 )->N();

            aJacobians( 1 )( tMasterNumDofDependencies )
            = - mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->testTraction( mNormal ) * mSlaveFI( 0 )->N()
              + mStabilizationParam( 0 )->val()( 0 ) * trans( mSlaveFI( 0 )->N() ) * mSlaveFI( 0 )->N();

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
                    += - trans( mMasterFI( 0 )->N() ) * mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal );
                       //+ mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump;

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mStabilizationParam( 1 )->val()( 0 ) * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal );
                }

                // if dependency of stabilization parameters on the dof type
                if ( mStabilizationParam( 0 )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tJump * mStabilizationParam( 0 )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += - trans( mSlaveFI( 0 )->N() ) * tJump * mStabilizationParam( 0 )->dSPdMasterDOF( tDofType );
                }

                if ( mStabilizationParam( 1 )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal ) * mStabilizationParam( 1 )->dSPdMasterDOF( tDofType )
                       + mMasterCM( 0 )->testTraction( mNormal ) * tJump * mStabilizationParam( 1 )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal ) * mStabilizationParam( 1 )->dSPdMasterDOF( tDofType );
                }

                if ( mStabilizationParam( 2 )->check_dof_dependency( tDofType, mtk::Master_Slave::MASTER ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mSlaveCM( 0 )->traction( mNormal ) * mStabilizationParam( 2 )->dSPdMasterDOF( tDofType );

                    aJacobians( 1 )( iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mSlaveCM( 0 )->traction( mNormal ) * mStabilizationParam( 2 )->dSPdMasterDOF( tDofType )
                       + mSlaveCM( 0 )->testTraction( mNormal ) * tJump * mStabilizationParam( 2 )->dSPdMasterDOF( tDofType );
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
                    += - trans( mMasterFI( 0 )->N() ) * mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->dTractiondDOF( tDofType, mNormal );
                       //+ mStabilizationParam( 2 )->val()( 0 ) * mSlaveCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tJump;
                }

                // if dependency of stabilization parameters on the dof type
                if ( mStabilizationParam( 0 )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tJump * mStabilizationParam( 0 )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mSlaveFI( 0 )->N() ) * tJump * mStabilizationParam( 0 )->dSPdSlaveDOF( tDofType );
                }

                if ( mStabilizationParam( 1 )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal ) * mStabilizationParam( 1 )->dSPdSlaveDOF( tDofType )
                       + mMasterCM( 0 )->testTraction( mNormal ) * tJump * mStabilizationParam( 1 )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += trans( mSlaveFI( 0 )->N() ) * mMasterCM( 0 )->traction( mNormal ) * mStabilizationParam( 1 )->dSPdSlaveDOF( tDofType );
                }

                if ( mStabilizationParam( 2 )->check_dof_dependency( tDofType, mtk::Master_Slave::SLAVE ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) * mSlaveCM( 0 )->traction( mNormal ) * mStabilizationParam( 2 )->dSPdSlaveDOF( tDofType );

                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                    +=   trans( mSlaveFI( 0 )->N() ) * mSlaveCM( 0 )->traction( mNormal ) * mStabilizationParam( 2 )->dSPdSlaveDOF( tDofType )
                       + mSlaveCM( 0 )->testTraction( mNormal ) * tJump * mStabilizationParam( 2 )->dSPdSlaveDOF( tDofType );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Interface::compute_jacobian_and_residual - This function does nothing.");
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
