
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Olsson_CLS_Bulk::IWG_Olsson_CLS_Bulk(       Field_Interpolator * aFieldInterpolator,
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

        void IWG_Olsson_CLS_Bulk::compute_residual( Matrix< DDRMat > & aResidual,
                                                    Matrix< DDRMat >   aFieldNormal )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            // evaluate the shape functions first space derivative
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // evaluate the shape function first time derivative
            Matrix< DDRMat > tBt = mFieldInterpolator->Bt();

            // evaluate the field value
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            //evaluate the field first space derivative
            Matrix< DDRMat > tPhiGradx = mFieldInterpolator->gradx( 1 );

            // evaluate the field first time derivative
            Matrix< DDRMat > tPhiGradt = mFieldInterpolator->gradt( 1 );

            //compute the residual
//            aResidual = trans( tN ) * tPhiGradt
//                      - trans( tBx ) * ( tPhi( 0 ) - mPhiLB ) * ( mPhiUB - tPhi( 0 ) ) * aFieldNormal
//                      + trans( tBx ) * mEpsilon * dot( tPhiGradx, aFieldNormal)        * aFieldNormal;

            aResidual = trans( tN ) * tPhiGradt
                      - trans( tBx ) * ( ( tPhi( 0 ) - mPhiLB ) * (mPhiUB - tPhi( 0 ) )
                                         - mEpsilon  * dot( tPhiGradx, aFieldNormal) ) * aFieldNormal;
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian( Matrix< DDRMat > & aJacobian,
                                                    Matrix< DDRMat >   aFieldNormal )
        {
            // evaluate the shape functions
            Matrix< DDRMat > tN = mFieldInterpolator->N();

            // evaluate the shape functions first space derivative
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();

            // evaluate the shape function first time derivative
            Matrix< DDRMat > tBt = mFieldInterpolator->Bt();

            // evaluate the field value
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            //compute the jacobian
            aJacobian = trans( tN )  * tBt
                      - trans( tBx ) * ( mPhiUB + mPhiLB - 2 * tPhi( 0 ) ) * aFieldNormal * tN
                      + trans( tBx ) * mEpsilon * ( aFieldNormal * trans( aFieldNormal ) * tBx );
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian_and_residual( Matrix< DDRMat > & aJacobian,
                                                                 Matrix< DDRMat > & aResidual,
                                                                 Matrix< DDRMat >   aFieldNormal )
        {
            // evaluate the shape functions and transpose
            Matrix< DDRMat > tN  = mFieldInterpolator->N();
            Matrix< DDRMat > tNt = trans( tN );

            // evaluate the shape functions first space derivative
            Matrix< DDRMat > tBx = mFieldInterpolator->Bx();
            Matrix< DDRMat > tBxt = trans( tBx );

            // evaluate the shape function first time derivative
            Matrix< DDRMat > tBt = mFieldInterpolator->Bt();

            // evaluate the field value
            Matrix< DDRMat > tPhi = mFieldInterpolator->val();

            //evaluate the field first space derivative
            Matrix< DDRMat > tPhiGradx = mFieldInterpolator->gradx( 1 );

            // evaluate the field first time derivative
            Matrix< DDRMat > tPhiGradt = mFieldInterpolator->gradt( 1 );

            //compute the residual
            aResidual = trans( tN )  * tPhiGradt
                      - trans( tBx ) * ( ( tPhi( 0 ) - mPhiLB ) * (mPhiUB - tPhi( 0 ) )
                                         - mEpsilon  * dot( tPhiGradx, aFieldNormal) ) * aFieldNormal;

            //compute the jacobian
            aJacobian = tNt  * tBt
                      - tBxt * ( ( mPhiUB + mPhiLB - 2 * tPhi( 0 ) ) * aFieldNormal * tN
                                 - mEpsilon * ( aFieldNormal * trans( aFieldNormal ) * tBx ) );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
