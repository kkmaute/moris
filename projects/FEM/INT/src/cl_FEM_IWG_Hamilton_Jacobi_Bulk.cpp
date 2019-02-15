
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Hamilton_Jacobi_Bulk::IWG_Hamilton_Jacobi_Bulk( Field_Interpolator * aFieldInterpolator )
        {
            // set the field interpolator
            mFieldInterpolator = aFieldInterpolator;
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_residual( Matrix< DDRMat > & aResidual,
                                                         Matrix< DDRMat >   aVelocityField )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            // evaluate the field first time derivative
            Matrix< DDRMat > tPhiGradt = mFieldInterpolator->gradt( 1 );

            // evaluate the field first space derivative
            Matrix< DDRMat > tPhiGradx = mFieldInterpolator->gradx( 1 );

           //compute the residual
            aResidual = trans( tN ) * ( tPhiGradt + trans( aVelocityField ) * tPhiGradx );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian( Matrix< DDRMat > & aJacobian,
                                                         Matrix< DDRMat >   aVelocityField )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            //evaluate the shape functions first derivative wrt t
            Matrix< DDRMat > tBt = mFieldInterpolator->Bt();

            // evaluate the shape functions first derivative wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // compute the jacobian
            aJacobian = trans( tN ) * ( tBt + trans( aVelocityField ) * tBx );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                                      Matrix< DDRMat > & aResidual,
                                                                      Matrix< DDRMat >   aVelocityField )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tNt  = mFieldInterpolator->N();
            tNt = trans( tNt );

            // evaluate the temporal gradient of the field
            Matrix< DDRMat > tphigradt = mFieldInterpolator->gradt( 1 );

            // evaluate the spatial gradient of the field
            Matrix< DDRMat > tphigradx = mFieldInterpolator->gradx( 1 );

            //evaluate the shape functions first derivative wrt t
            Matrix< DDRMat > tBt = mFieldInterpolator->Bt();

            // evaluate the shape functions first derivative wrt x
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // transpose the velocity field
            aVelocityField = trans( aVelocityField );

            //compute the residual
            aResidual = tNt * ( tphigradt + aVelocityField * tphigradx );

            // compute the jacobian
            aJacobian = tNt * ( tBt +  aVelocityField * tBx );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
