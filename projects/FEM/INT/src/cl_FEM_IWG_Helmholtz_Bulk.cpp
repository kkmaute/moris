
#include "cl_FEM_IWG_Helmholtz_Bulk.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------
        IWG_Helmholtz_Bulk::IWG_Helmholtz_Bulk()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

//            // set the residual dof type
//            mResidualDofType = { MSI::Dof_Type::VX };
//
//            // set the active dof type
//            mMasterDofTypes = {{ MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            //FIXME set unfiltered velocity values at nodes
            Matrix< DDRMat > tVHat  = mNodalWeakBCs;

            // set field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->gradx( 1 )
                           + trans( vN->N() ) * ( vN->val() - vN->N() * tVHat );
        }

//------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // set field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian
            aJacobians( 0 )( 0 ) = mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->dnNdxn( 1 )
                                 + trans( vN->N() ) * vN->N();
        }

//------------------------------------------------------------------------------
        void IWG_Helmholtz_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            //FIXME set unfiltered velocity values at nodes
            Matrix< DDRMat > tVHat  = mNodalWeakBCs;

            // set field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            // set residual size
            this->set_residual( aResidual );

            // compute the residual
            aResidual( 0 ) = mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->gradx( 1 )
                           + trans( vN->N() ) * ( vN->val() - vN->N() * tVHat );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the residual
            aJacobians( 0 )( 0 ) = mFilterParam * trans( vN->dnNdxn( 1 ) ) * vN->dnNdxn( 1 )
                                 + trans( vN->N() ) * vN->N();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
