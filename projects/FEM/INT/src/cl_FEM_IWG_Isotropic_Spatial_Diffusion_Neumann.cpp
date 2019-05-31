
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

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_residual
            ( Matrix< DDRMat >                   & aResidual,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the normal flux
            Matrix < DDRMat > tNormalFlux = tTemp->N() * mNodalWeakBCs;

            // compute the residual r_T
            aResidual = - trans( tTemp->N() ) * tNormalFlux;
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian
            ( moris::Cell< Matrix< DDRMat > >    & aJacobians,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian_and_residual
            ( moris::Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >                   & aResidual,
              moris::Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* tTemp = aFieldInterpolators( 0 );

            // compute the residual r_T
            aResidual = - trans( tTemp->N() ) * tTemp->N() * mNodalWeakBCs;

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian j_T_T
            uint tNumOfBases = tTemp->get_number_of_space_time_bases();
            aJacobians( 0 ).set_size( tNumOfBases, tNumOfBases, 0.0 );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
