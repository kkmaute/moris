
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

            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute flux
            Matrix< DDRMat > tStress;
            mMasterCM( 0 )->eval_stress( tStress );

            // compute the residual
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * tStress - trans( mMasterFI( 0 )->N() ) * tQ;
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

            // loop over the dof dependencies
            uint tNumDofTypes = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // compute flux derivative
                Matrix< DDRMat > tdStressdDOF;
                mMasterCM( 0 )->eval_dStressdDOF( mMasterGlobalDofTypes( iDOF ), tdStressdDOF );

                // compute the jacobian
                aJacobians( 0 )( iDOF ) = trans( mMasterFI( 0 )->Bx() ) * tdStressdDOF;
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            //FIXME heat load enforced
            Matrix< DDRMat > tQ( 1, 1, 0.0 );

            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the residual size
            this->set_residual( aResidual );

            // compute flux
            Matrix< DDRMat > tStress;
            mMasterCM( 0 )->eval_stress( tStress );

            // compute the residual
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * tStress - trans( mMasterFI( 0 )->N() ) * tQ;

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // loop over dof dependencies
            uint tNumDofTypes = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofTypes; iDOF++ )
            {
                // compute flux derivative
                Matrix< DDRMat > tdStressdDOF;
                mMasterCM( 0 )->eval_dStressdDOF( mMasterGlobalDofTypes( iDOF ), tdStressdDOF );

                // compute jacobian
                aJacobians( 0 )( iDOF ) = trans( mMasterFI( 0 )->Bx() ) * tdStressdDOF;
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
