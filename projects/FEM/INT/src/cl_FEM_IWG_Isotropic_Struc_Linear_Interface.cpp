
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Interface.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
    IWG_Isotropic_Struc_Linear_Interface::IWG_Isotropic_Struc_Linear_Interface()
        {
            // FIXME set penalty parameter
            mGammaInterface = 1.0;

            // FIXME set weight master parameter
            mMasterWeight = 0.5;

            // FIXME set weight slave parameter
            mSlaveWeight = 0.5;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
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
            Matrix< DDRMat > tMasterFlux;
            mMasterCM( 0 )->eval_flux( tMasterFlux );

            // evaluate slave flux
            Matrix< DDRMat > tSlaveFlux;
            mSlaveCM( 0 )->eval_flux( tSlaveFlux );

            // evaluate average flux
            Matrix< DDRMat > tFlux = mMasterWeight * tMasterFlux + mSlaveWeight * tSlaveFlux;

            // evaluate master conductivity matrix
            Matrix< DDRMat > tMasterK;
            mMasterCM( 0 )->eval_const( tMasterK );

            // evaluate slave conductivity matrix
            Matrix< DDRMat > tSlaveK;
            mSlaveCM( 0 )->eval_const( tSlaveK );

            // evaluate temperature jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mSlaveFI( 0 )->val();

            Matrix< DDRMat > tMasterTestStrain;
            mMasterCM( 0 )->eval_test_strain( tMasterTestStrain );

            Matrix< DDRMat > tSlaveTestStrain;
            mSlaveCM( 0 )->eval_test_strain( tSlaveTestStrain );

            Matrix< DDRMat > tNormal(2 ,3, 0.0 );
            tNormal( 0, 0 ) = mNormal( 0,0 );
            tNormal( 0, 2 ) = mNormal( 1,0 );
            tNormal( 1, 1 ) = mNormal( 1,0 );
            tNormal( 1, 2 ) = mNormal( 0,0 );

            Matrix< DDRMat >tShapeFuncMaster2D( 2, 8, 0.0 );
            tShapeFuncMaster2D( 0, 0 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFuncMaster2D( 0, 1 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFuncMaster2D( 0, 2 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFuncMaster2D( 0, 3 ) = mMasterFI( 0 )->N()( 0, 3 );
            tShapeFuncMaster2D( 1, 4 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFuncMaster2D( 1, 5 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFuncMaster2D( 1, 6 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFuncMaster2D( 1, 7 ) = mMasterFI( 0 )->N()( 0, 3 );

            Matrix< DDRMat >tShapeFuncSlave2D( 2, 8, 0.0 );
            tShapeFuncSlave2D( 0, 0 ) = mSlaveFI( 0 )->N()( 0, 0 );
            tShapeFuncSlave2D( 0, 1 ) = mSlaveFI( 0 )->N()( 0, 1 );
            tShapeFuncSlave2D( 0, 2 ) = mSlaveFI( 0 )->N()( 0, 2 );
            tShapeFuncSlave2D( 0, 3 ) = mSlaveFI( 0 )->N()( 0, 3 );
            tShapeFuncSlave2D( 1, 4 ) = mSlaveFI( 0 )->N()( 0, 0 );
            tShapeFuncSlave2D( 1, 5 ) = mSlaveFI( 0 )->N()( 0, 1 );
            tShapeFuncSlave2D( 1, 6 ) = mSlaveFI( 0 )->N()( 0, 2 );
            tShapeFuncSlave2D( 1, 7 ) = mSlaveFI( 0 )->N()( 0, 3 );

//            print( tMasterTestStrain ,"tMasterTestStrain");
//            print( tSlaveTestStrain  ,"tSlaveTestStrain");

            // compute master residual
            aResidual( 0 ) = - trans( tShapeFuncMaster2D ) * tNormal * tFlux
                             + mMasterWeight * trans( tMasterTestStrain ) * tMasterK * trans( tNormal ) * trans(tJump)
                             + mGammaInterface * trans( tShapeFuncMaster2D ) * trans( tJump );

            // compute slave residual
            aResidual( 1 ) =   trans( tShapeFuncSlave2D ) * tNormal * tFlux
                             + mSlaveWeight * trans( tSlaveTestStrain ) * tSlaveK * trans( tNormal ) * trans(tJump)
                             - mGammaInterface * trans( tShapeFuncSlave2D ) * trans( tJump );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Interface::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
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

            Matrix< DDRMat > tMasterTestStrain;
            mMasterCM( 0 )->eval_test_strain( tMasterTestStrain );

            Matrix< DDRMat > tSlaveTestStrain;
            mSlaveCM( 0 )->eval_test_strain( tSlaveTestStrain );

            Matrix< DDRMat > tNormal(2 ,3, 0.0 );
            tNormal( 0, 0 ) = mNormal( 0,0 );
            tNormal( 0, 2 ) = mNormal( 1,0 );
            tNormal( 1, 1 ) = mNormal( 1,0 );
            tNormal( 1, 2 ) = mNormal( 0,0 );

            Matrix< DDRMat >tShapeFuncMaster2D( 2, 8, 0.0 );
            tShapeFuncMaster2D( 0, 0 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFuncMaster2D( 0, 1 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFuncMaster2D( 0, 2 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFuncMaster2D( 0, 3 ) = mMasterFI( 0 )->N()( 0, 3 );
            tShapeFuncMaster2D( 1, 4 ) = mMasterFI( 0 )->N()( 0, 0 );
            tShapeFuncMaster2D( 1, 5 ) = mMasterFI( 0 )->N()( 0, 1 );
            tShapeFuncMaster2D( 1, 6 ) = mMasterFI( 0 )->N()( 0, 2 );
            tShapeFuncMaster2D( 1, 7 ) = mMasterFI( 0 )->N()( 0, 3 );

            Matrix< DDRMat >tShapeFuncSlave2D( 2, 8, 0.0 );
            tShapeFuncSlave2D( 0, 0 ) = mSlaveFI( 0 )->N()( 0, 0 );
            tShapeFuncSlave2D( 0, 1 ) = mSlaveFI( 0 )->N()( 0, 1 );
            tShapeFuncSlave2D( 0, 2 ) = mSlaveFI( 0 )->N()( 0, 2 );
            tShapeFuncSlave2D( 0, 3 ) = mSlaveFI( 0 )->N()( 0, 3 );
            tShapeFuncSlave2D( 1, 4 ) = mSlaveFI( 0 )->N()( 0, 0 );
            tShapeFuncSlave2D( 1, 5 ) = mSlaveFI( 0 )->N()( 0, 1 );
            tShapeFuncSlave2D( 1, 6 ) = mSlaveFI( 0 )->N()( 0, 2 );
            tShapeFuncSlave2D( 1, 7 ) = mSlaveFI( 0 )->N()( 0, 3 );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterWeight * trans( tMasterTestStrain ) * tMasterK * trans( tNormal ) * tShapeFuncMaster2D
                                   + mGammaInterface * trans( tShapeFuncMaster2D ) * tShapeFuncMaster2D;

            aJacobians( 0 )( 1 ) = - mMasterWeight * trans( tMasterTestStrain ) * tMasterK * trans( tNormal ) * tShapeFuncSlave2D
                                   - mGammaInterface * trans( tShapeFuncMaster2D ) * tShapeFuncSlave2D;

            aJacobians( 1 )( 0 ) = + mSlaveWeight * trans( tSlaveTestStrain ) * tSlaveK * trans( tNormal ) * tShapeFuncMaster2D
                                   - mGammaInterface * trans( tShapeFuncSlave2D ) * tShapeFuncMaster2D;

            aJacobians( 1 )( 1 ) = - mSlaveWeight * trans( tSlaveTestStrain ) * tSlaveK * trans( tNormal ) * tShapeFuncSlave2D
                                   + mGammaInterface * trans( tShapeFuncSlave2D ) * tShapeFuncSlave2D;

            // compute the jacobian for indirect dof dependencies through master constitutive models
            uint tMasterNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tMasterNumDofDependencies; iDOF++ )
            {
                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // evaluate master flux
                    Matrix< DDRMat > tMasterdFlux;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tMasterdFlux );

                    // evaluate master contitutive matrix
                    Matrix< DDRMat > tMasterdConst;
                    mMasterCM( 0 )->eval_dConstdDOF( mMasterGlobalDofTypes( iDOF ), tMasterdConst );

                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                        += - mMasterWeight * trans( tShapeFuncMaster2D ) *  tNormal * tMasterdFlux;
                           //+ mMasterWeight * trans( tMasterTestStrain ) * tMasterdConst * trans( tNormal) * tJump;
                    aJacobians( 1 )( iDOF ).matrix_data()
                        += mMasterWeight * trans( tShapeFuncSlave2D ) * tNormal * tMasterdFlux;
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
                    Matrix< DDRMat > tSlavedFlux;
                    mSlaveCM( 0 )->eval_dFluxdDOF( mSlaveGlobalDofTypes( iDOF ), tSlavedFlux );

                    // evaluate slave contitutive matrix
                    Matrix< DDRMat > tSlavedConst;
                    mSlaveCM( 0 )->eval_dConstdDOF( mSlaveGlobalDofTypes( iDOF ), tSlavedConst );

                    // add contribution to jacobian
                    aJacobians( 0 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                        += - mSlaveWeight * trans( tShapeFuncMaster2D ) * tNormal * tSlavedFlux;
                    aJacobians( 1 )( tMasterNumDofDependencies + iDOF ).matrix_data()
                        += mSlaveWeight * trans( tShapeFuncSlave2D ) *  tNormal * tSlavedFlux;
                         //+ mSlaveWeight * trans( tSlaveTestStrain ) * tSlavedConst * trans( tNormal ) * tJump;
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
