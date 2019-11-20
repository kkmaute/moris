
#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // FIXME check if selection matrix set, if not set to identity
            Matrix< DDRMat > tM;
            if( mMasterProp.size() > 1 )
            {
                tM = mMasterProp( 1 )->val();
            }
            else
            {
                uint tSpaceDim = mMasterFI( 0 )->get_dof_type().size();
                eye( tSpaceDim, tSpaceDim, tM );
            }

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * tM * mMasterCM( 0 )->traction( mNormal )
                             + mMasterCM( 0 )->testTraction( mNormal ) * tM * tJump
                             + mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * tJump;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

             // FIXME check if selection matrix was set, if not set to identity
             Matrix< DDRMat > tM;
             if( mMasterProp.size() > 1 )
             {
                 tM = mMasterProp( 1 )->val();
             }
             else
             {
                 uint tSpaceDim = mMasterFI( 0 )->get_dof_type().size();
                 eye( tSpaceDim, tSpaceDim, tM );
             }

            // compute jump
            Matrix< DDRMat > tJump = mMasterFI( 0 )->val() - mMasterProp( 0 )->val();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            aJacobians( 0 )( 0 ) = mMasterCM( 0 )->testTraction( mNormal ) * tM * mMasterFI( 0 )->N()
                                 + mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * mMasterFI( 0 )->N();

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
                    += -1.0 * mMasterCM( 0 )->testTraction( mNormal ) * tM * mMasterProp( 0 )->dPropdDOF( tDofType )
                       - mStabilizationParam( 0 )->val()( 0 ) * trans( mMasterFI( 0 )->N() ) * tM * mMasterProp( 0 )->dPropdDOF( tDofType );
                }

                // if dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += - trans( mMasterFI( 0 )->N() ) *  tM * mMasterCM( 0 )->dTractiondDOF( tDofType, mNormal ) ;
                    // + mMasterCM( 0 )->dTestTractiondDOF( tDofType, mNormal ) * tM * tJump;
                }

                // if dependency on the dof type
                if ( mStabilizationParam( 0 )->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterFI( 0 )->N() ) * tM * tJump * mStabilizationParam( 0 )->dSPdMasterDOF( tDofType );
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
