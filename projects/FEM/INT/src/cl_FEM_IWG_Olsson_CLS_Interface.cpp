
#include "cl_FEM_IWG_Olsson_CLS_Interface.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Olsson_CLS_Interface::IWG_Olsson_CLS_Interface(       Field_Interpolator * aFieldInterpolator,
                                                            const real                 aFieldUpperBound,
                                                            const real                 aFieldLowerBound,
                                                            const real                 aMethodParameter )
        {
            // set the field interpolator
            mFieldInterpolator = aFieldInterpolator;

            //set field upper and lower bound
            mPhiUB = aFieldUpperBound;
            mPhiLB = aFieldLowerBound;

            // set Olsson CLS epsilon parameter
            mEpsilon = aMethodParameter;

        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Interface::compute_residual( Matrix< DDRMat > & aResidual,
                                                         Matrix< DDRMat >   aFieldNormal,
                                                         Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tN = mFieldInterpolator->N();
            Matrix< DDRMat > tNt = trans( tN );

            // evaluate the field
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            // evaluate the field first space derivative
            Matrix< DDRMat > tPhiGradx = mFieldInterpolator->gradx( 1 );

            //compute the residual
            aResidual = tNt * ( ( tPhi( 0 ) - mPhiLB ) * ( mPhiUB - tPhi( 0 ) )
                                - mEpsilon * dot( tPhiGradx, aFieldNormal ) )
                            * trans( aFieldNormal ) * aInterfaceNormal;
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Interface::compute_jacobian( Matrix< DDRMat > & aJacobian,
                                                         Matrix< DDRMat >   aFieldNormal,
                                                         Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tN = mFieldInterpolator->N();
            Matrix< DDRMat > tNt = trans( tN );

            // evaluate the shape function first space derivative
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // evaluate the field
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            // compute the jacobian
            aJacobian = tNt * ( ( mPhiLB + mPhiUB - 2 * tPhi ) * tN
                                - mEpsilon * trans( aFieldNormal) * tBx ) * dot( aFieldNormal, aInterfaceNormal ) ;

        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Interface::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                                      Matrix< DDRMat > & aResidual,
                                                                      Matrix< DDRMat >   aFieldNormal,
                                                                      Matrix< DDRMat >   aInterfaceNormal )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tN = mFieldInterpolator->N();
            Matrix< DDRMat > tNt = trans( tN );

            // evaluate the shape function first space derivative
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // evaluate the field
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            // evaluate the field first space derivative
            Matrix< DDRMat > tPhiGradx = mFieldInterpolator->gradx( 1 );

            //compute the residual
            aResidual = tNt * ( ( tPhi( 0 ) - mPhiLB ) * ( mPhiUB - tPhi( 0 ) )
                                - mEpsilon * dot( tPhiGradx, aFieldNormal ) )
                            * trans( aFieldNormal ) * aInterfaceNormal;

            // compute the jacobian
            aJacobian = tNt * ( ( mPhiLB + mPhiUB - 2 * tPhi ) * tN
                                - mEpsilon * trans( aFieldNormal) * tBx ) * dot( aFieldNormal, aInterfaceNormal );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
