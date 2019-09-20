
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val();
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct IWG dof dependencies
            uint tNumOfBases = mMasterFI( 0 )->get_number_of_space_time_bases();
            aJacobians( 0 )( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // if dependency in the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( 0 ) );
                }
            }
        }

//------------------------------------------------------------------------------
        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                                     moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // set residual size
            this->set_residual( aResidual );

            // compute the residual r_T
            aResidual( 0 ) = - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->val();

            // set jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian for direct IWG dof dependencies
            uint tNumOfBases = mMasterFI( 0 )->get_number_of_space_time_bases();
            aJacobians( 0 )( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );

            // compute the jacobian for indirect IWG dof dependencies through properties
            for( uint iDOF = 0; iDOF < mMasterGlobalDofTypes.size(); iDOF++ )
            {
                // if dependency in the dof type
                if ( mMasterProp( 0 )->check_dof_dependency( mMasterGlobalDofTypes( iDOF ) ) )
                {
                    // add contribution to jacobian
                    aJacobians( 0 )( iDOF ).matrix_data() += - trans( mMasterFI( 0 )->N() ) * mMasterProp( 0 )->dPropdDOF( mMasterGlobalDofTypes( 0 ) );
                }
            }
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
