
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
            //FIXME: velocity scalar field not UX
            mResidualDofType = MSI::Dof_Type::UX;

            // set the active dof type
            //FIXME: velocity scalar field not UX
            mActiveDofTypes = {{ MSI::Dof_Type::UX }};
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_residual( Matrix< DDRMat >            & aResidual,
                                                        Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set the field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // compute the residual
            aResidual = - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                        Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set the field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // compute the jacobian
            aJacobians( 0 ) = - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->Bx();
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Interface::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                     Matrix< DDRMat >            & aResidual,
                                                                     Cell< Field_Interpolator* > & aFieldInterpolators)
        {
            // set the field interpolator
            Field_Interpolator* vN = aFieldInterpolators( 0 );

            //FIXME set the interface normal
            Matrix< DDRMat > aInterfaceNormal( vN->gradx( 1 ).n_cols() , 1, 1.0 );

            // compute the residual
            aResidual = - mFilterParam * trans( vN->N() ) * trans( vN->gradx( 1 ) ) * aInterfaceNormal;

            // compute the residual
            aJacobians( 0 ) = - mFilterParam * trans( vN->N() ) * trans( aInterfaceNormal ) * vN->Bx();
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
