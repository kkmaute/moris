
#include "cl_FEM_IWG_Helmholtz_Interface.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//        IWG_Helmholtz_Interface::IWG_Helmholtz_Interface( const real aFilterParam )
        IWG_Helmholtz_Interface::IWG_Helmholtz_Interface()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::VX };

            // set the active dof type
            mActiveDofTypes = {{ MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_residual( Matrix< DDRMat >                   & aResidual )
        {
            // set the field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // compute the residual
            aResidual = - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians )
        {
            // set the field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian
            aJacobians( 0 ) = - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->Bx();
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                                     Matrix< DDRMat >                   & aResidual )
        {
            // set the field interpolator
            Field_Interpolator* vN = mMasterFI( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // compute the residual
            aResidual = - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal;

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the residual
            aJacobians( 0 ) = - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->Bx();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
