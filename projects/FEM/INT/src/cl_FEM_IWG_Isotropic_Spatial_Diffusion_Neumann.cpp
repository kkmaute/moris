
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
        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_residual
            ( Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {

        }

//------------------------------------------------------------------------------

        void IWG_Isotropic_Spatial_Diffusion_Neumann::compute_jacobian_and_residual
            ( Cell< Matrix< DDRMat > >    & aJacobians,
              Matrix< DDRMat >            & aResidual,
              Cell< Field_Interpolator* > & aFieldInterpolators )
        {

        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
