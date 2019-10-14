
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
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ).matrix_data() += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->flux();

            // if body load
            if ( mMasterPropTypes( 0 ) == fem::Property_Type::TEMP_LOAD )
            {
                aResidual( 0 ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val()( 0 );
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
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
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if we have a body load
                if ( mMasterPropTypes( 0 ) == fem::Property_Type::TEMP_LOAD )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->dFluxdDOF( tDofType );
                    // fixme add derivative of the test strain
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();
            this->check_properties();
            this->check_constitutive_models();

            // set the residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ).matrix_data() += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->flux();

            // compute residual due to body load
            if ( mMasterProp( 0 )->get_property_type() == fem::Property_Type::TEMP_LOAD )
            {
                aResidual( 0 ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val()( 0 );
            }

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for indirect dof dependencies through properties and constitutive model
            uint tNumDofDependencies = mMasterGlobalDofTypes.size();
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                Cell< MSI::Dof_Type > tDofType = mMasterGlobalDofTypes( iDOF );

                // if we have a body load
                if ( mMasterProp( 0 )->get_property_type() == fem::Property_Type::TEMP_LOAD )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( 0 )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( tDofType );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( 0 )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterCM( 0 )->testStrain() ) * mMasterCM( 0 )->dFluxdDOF( tDofType );
                    // fixme add test strain derivative
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
