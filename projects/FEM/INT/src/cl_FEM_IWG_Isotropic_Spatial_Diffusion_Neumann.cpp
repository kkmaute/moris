
#include "cl_FEM_IWG_Isotropic_Spatial_Diffusion_Neumann.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Isotropic_Spatial_Diffusion_Neumann::IWG_Isotropic_Spatial_Diffusion_Neumann()
        {

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::TEMP };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::TEMP } };

            // set the active mp type
            mActivePropertyTypes = { fem::Property_Type::TEMP_NEUMANN };
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_residual( Matrix< DDRMat > & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            //this->check_properties();

            // fixme compute the normal flux
            Matrix < DDRMat > tNormalFlux = mMasterFI( 0 )->N() * mNodalWeakBCs;
            //Matrix < DDRMat > tNormalFlux;
            //mMasterProp( 0 )->val( tNormalFlux );

            // compute the residual r_T
            aResidual = - trans( mMasterFI( 0 )->N() ) * tNormalFlux;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian( moris::Cell< Matrix< DDRMat > > & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = mMasterFI( 0 )->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > > & aJacobians,
                                                                                     Matrix< DDRMat >                & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators();

            // check master properties
            this->check_properties();

            // fixme compute the normal flux
            Matrix < DDRMat > tNormalFlux = mMasterFI( 0 )->N() * mNodalWeakBCs;
            //Matrix < DDRMat > tNormalFlux;
            //mMasterProp( 0 )->val( tNormalFlux );

            // compute the residual r_T
            aResidual = - trans( mMasterFI( 0 )->N() ) * tNormalFlux;

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = mMasterFI( 0 )->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
