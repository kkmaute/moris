
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk_Test.hpp"
#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        IWG_Hamilton_Jacobi_Bulk_Test::IWG_Hamilton_Jacobi_Bulk_Test()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS1 };

            // set the active dof type
            mMasterDofTypes = {{ MSI::Dof_Type::LS1 }};

        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );

            // velocity field value
            Matrix< DDRMat > aVN( 1, 3, 1.0 );

            // set residual size
            this->set_residual( aResidual );

           //compute the residual
           aResidual( 0 ) = trans( phi->N() ) * ( phi->gradt( 1 ) + aVN * phi->gradx( 1 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );

            // velocity field value
            Matrix< DDRMat > aVN( 1, 3, 1.0 );

            // set the jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian Jphiphi
            aJacobians( 0 )( 0 ) = trans( phi->N() ) * ( phi->Bt() + aVN * phi->Bx() );

        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk_Test::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                           moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );

            // velocity field value
            Matrix< DDRMat > aVN( 1, 3, 1.0 );

            // set the jacobian size
            this->set_residual( aResidual );

            //compute the residual
            aResidual( 0 ) = trans( phi->N() ) * ( phi->gradt( 1 ) + aVN * phi->gradx( 1 ) );

            // set the jacobian size
            this->set_jacobian( aJacobians );


            // compute the jacobian Jphiphi
            aJacobians( 0 )( 0 ) = trans( phi->N() ) * ( phi->Bt() + aVN * phi->Bx() );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
