
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Helmholtz_Bulk::IWG_Helmholtz_Bulk(       Field_Interpolator * aFieldInterpolator,
                                                const real                 aFilterParam )
        {
            // set the field interpolator
            mFieldInterpolator = aFieldInterpolator;

            //set the Helmholtz filter parameter
            mFilterParam       = aFilterParam;

        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_residual( Matrix< DDRMat > & aResidual,
                                                   real               aVHat)
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the shape function first derivatives wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            //evaluate the spatial gradient of the field
            Matrix< DDRMat > tvgradx = mFieldInterpolator->gradx(1);

            //evaluate the field
            Matrix< DDRMat > tv = mFieldInterpolator->val();

            // compute the residual
            aResidual = mFilterParam * trans( tBx ) * tvgradx
                      + trans( tN ) * tv
                      - trans( tN ) * aVHat;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_jacobian( Matrix< DDRMat > & aJacobian )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the shape function first derivatives wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // compute the jacobian
            aJacobian = mFilterParam * trans( tBx ) * tBx
                      + trans( tN ) * tN;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                                Matrix< DDRMat > & aResidual,
                                                                real               aVHat )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the shape function first derivatives wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            //evaluate the spatial gradient of the field
            Matrix< DDRMat > tvgradx = mFieldInterpolator->gradx(1);

            //evaluate the field
            Matrix< DDRMat > tv = mFieldInterpolator->val();

            // compute the residual
            aResidual = mFilterParam * trans( tBx ) * tvgradx
                      + trans( tN ) * tv
                      - trans( tN ) * aVHat;

            // compute the residual
            aJacobian = mFilterParam * trans( tBx ) * tBx
                      + trans( tN ) * tN;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
