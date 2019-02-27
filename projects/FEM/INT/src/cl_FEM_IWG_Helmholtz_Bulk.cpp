
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"
#include "fn_trans.hpp"
//#include "op_times.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//        IWG_Helmholtz_Bulk::IWG_Helmholtz_Bulk( const real aFilterParam )
        IWG_Helmholtz_Bulk::IWG_Helmholtz_Bulk()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            // set the residual dof type
            mResidualDofType = MSI::Dof_Type::VX;

            // set the active dof type
            mActiveDofTypes = {{ MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_residual( Matrix< DDRMat >            & aResidual,
                                                   Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            //FIXME set unfiltered velocity value
            real aVHat  = 1;

            // set field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            // compute the residual
            aResidual = mFilterParam * trans( vN->Bx() ) * vN->gradx( 1 )
                      + trans( vN->N() ) * ( vN->val() - aVHat );
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                   Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            // compute the jacobian
            aJacobians( 0 ) = mFilterParam * trans( vN->Bx() ) * vN->Bx()
                            + trans( vN->N() ) * vN->N();
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                Matrix< DDRMat >            & aResidual,
                                                                Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            //FIXME set unfiltered velocity value
            real aVHat  = 1;

            // set field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            // compute the residual
            aResidual = mFilterParam * trans( vN->Bx() ) * vN->gradx( 1 )
                      + trans( vN->N() ) * ( vN->val() - aVHat );

            // compute the residual
            aJacobians( 0 ) = mFilterParam * trans( vN->Bx() ) * vN->Bx()
                            + trans( vN->N() ) * vN->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
