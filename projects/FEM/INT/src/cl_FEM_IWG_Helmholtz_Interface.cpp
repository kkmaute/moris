
#include "cl_FEM_IWG_Helmholtz_Interface.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Helmholtz_Interface::IWG_Helmholtz_Interface(       Field_Interpolator * aFieldInterpolator,
                                                          const real                 aFilterParam )
        {
            // set the field interpolator
            mFieldInterpolator = aFieldInterpolator;

            //set the Helmholtz filter parameter
            mFilterParam       = aFilterParam;

        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_residual( Matrix< DDRMat > & aResidual,
                                                        Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the spatial gradient of the field
            Matrix< DDRMat > tVgradx = mFieldInterpolator->gradx(1);

            // compute the residual
            aResidual = - mFilterParam * trans( tN ) * trans( tVgradx ) * aInterfaceNormal;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian( Matrix< DDRMat > & aJacobian,
                                                        Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the shape function first derivatives wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // compute the jacobian
            aJacobian = - mFilterParam * trans( tN ) * trans( aInterfaceNormal ) * tBx;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                                     Matrix< DDRMat > & aResidual,
                                                                     Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tNt = mFieldInterpolator->N();
            tNt = trans( tNt );

            //evaluate the shape function first derivatives wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            //evaluate the spatial gradient of the field
            Matrix< DDRMat > tVgradx = mFieldInterpolator->gradx(1);

            // compute the residual
            aResidual = - mFilterParam * tNt * trans( tVgradx ) * aInterfaceNormal;

            // compute the residual
            aJacobian = - mFilterParam * tNt * trans( aInterfaceNormal ) * tBx;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
