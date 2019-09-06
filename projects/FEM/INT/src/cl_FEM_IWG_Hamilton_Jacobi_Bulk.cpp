
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk.hpp"

#include "fn_trans.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        IWG_Hamilton_Jacobi_Bulk::IWG_Hamilton_Jacobi_Bulk()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS1 };

            // set the active dof type
            mMasterDofTypes = {{ MSI::Dof_Type::LS1 },
                               { MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_residual( moris::Cell< Matrix< DDRMat > > & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

            // set residual size
            this->set_residual( aResidual );

           //compute the residual
           aResidual( 0 ) = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * phi->gradx( 1 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

            // set jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian Jphiphi
            aJacobians( 0 )( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val() * phi->Bx() );

            // compute the jacobian JphivN
            uint tvNNumOfDofs = vN->get_number_of_fields()*vN->get_number_of_space_time_bases();
            aJacobians( 0 )( 1 ) = trans( phi->N() ) * reshape( trans( vN->N() ) * trans( phi->gradx( 1 ) ), 1, tvNNumOfDofs );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian_and_residual( moris::Cell< moris::Cell< Matrix< DDRMat > > > & aJacobians,
                                                                      moris::Cell< Matrix< DDRMat > >                & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

            // set residual size
            this->set_residual( aResidual );

            //compute the residual
            aResidual( 0 ) = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * phi->gradx( 1 ) );

            // set jacobian size
            this->set_jacobian( aJacobians );

            // compute the jacobian Jphiphi
            aJacobians( 0 )( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val() * phi->Bx() );

            // compute the jacobian JphivN
            uint tvNNumOfDofs = vN->get_number_of_fields()*vN->get_number_of_space_time_bases();
            aJacobians( 0 )( 1 ) = trans( phi->N() ) * reshape( trans( vN->N() ) * trans( phi->gradx( 1 ) ), 1, tvNNumOfDofs );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
