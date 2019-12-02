
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------
        IWG_Isotropic_Spatial_Diffusion_Bulk::IWG_Isotropic_Spatial_Diffusion_Bulk()
        {
            // set size for the property pointer cell
            mMasterProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

            // populate the property map
            mPropertyMap[ "Load" ] = IWG_Property_Type::LOAD;

            // set size for the constitutive model pointer cell
            mMasterCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

            // populate the constitutive map
            mConstitutiveMap[ "DiffLinIso" ] = IWG_Constitutive_Type::DIFF_LIN_ISO;
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ).matrix_data() += trans( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->testStrain() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO) )->flux();

            // if body load
            if ( mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) ) != nullptr )
            {
                aResidual( 0 ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD) )->val()( 0 );
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators, properties and constitutive models
            this->check_dof_field_interpolators();
            this->check_dv_field_interpolators();

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
                if ( mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) ) != nullptr )
                {
                    // if property has dependency on the dof type
                    if ( mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) )->check_dof_dependency( tDofType ) )
                    {
                        // compute the jacobian
                        aJacobians( 0 )( iDOF ).matrix_data()
                        += - trans( mMasterFI( 0 )->N() ) * mMasterProp( static_cast< uint >( IWG_Property_Type::LOAD ) )->dPropdDOF( tDofType );
                    }
                }

                // if constitutive model has dependency on the dof type
                if ( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->check_dof_dependency( tDofType ) )
                {
                    // compute the jacobian
                    aJacobians( 0 )( iDOF ).matrix_data()
                    += trans( mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->testStrain() ) * mMasterCM( static_cast< uint >( IWG_Constitutive_Type::DIFF_LIN_ISO ) )->dFluxdDOF( tDofType );
                    // fixme add derivative of the test strain
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                  moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
