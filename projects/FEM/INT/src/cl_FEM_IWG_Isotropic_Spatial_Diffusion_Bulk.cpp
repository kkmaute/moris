
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
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute flux
            Matrix< DDRMat > tFlux;
            mMasterCM( 0 )->eval_flux( tFlux );

            // compute the residual
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * tFlux
                           - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val()( 0 );
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct dof dependencies
            // Here no direct dependencies

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // if property has dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute flux derivative
                    Matrix< DDRMat > tdFluxdDOF;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdFluxdDOF );

                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data() += trans( mMasterFI( 0 )->dnNdxn( 1 ) ) * tdFluxdDOF;
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the residual size
            this->set_residual( aResidual );

            // compute flux
            Matrix< DDRMat > tFlux;
            mMasterCM( 0 )->eval_flux( tFlux );

            // compute the residual
            aResidual( 0 ) = trans( mMasterFI( 0 )->Bx() ) * tFlux
                           - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val();

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // if property has dependency on the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( iDOF ) );
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // compute flux derivative
                    Matrix< DDRMat > tdFluxdDOF;
                    mMasterCM( 0 )->eval_dFluxdDOF( mMasterGlobalDofTypes( iDOF ), tdFluxdDOF );

                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data() += trans( mMasterFI( 0 )->Bx() ) * tdFluxdDOF;
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
