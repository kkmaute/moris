
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
            mActiveDofTypes = {{ MSI::Dof_Type::LS1 },
                               { MSI::Dof_Type::VX }};
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_residual( Matrix< DDRMat >            & aResidual,
                                                         Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolators
            Field_Interpolator* phi = aFieldInterpolators( 0 );
            Field_Interpolator* vN  = aFieldInterpolators( 1 );

           //compute the residual
           aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * phi->gradx( 1 ) );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                         Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolators
            Field_Interpolator* phi = aFieldInterpolators( 0 );
            Field_Interpolator* vN  = aFieldInterpolators( 1 );

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute the jacobian Jphiphi
            aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val() * phi->Bx() );

            // compute the jacobian JphivN
            uint tvNNumOfDofs = vN->get_number_of_fields()*vN->get_number_of_space_time_bases();
            aJacobians( 1 ) = trans( phi->N() ) * reshape( trans( vN->N() ) * trans( phi->gradx( 1 ) ), 1, tvNNumOfDofs );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                      Matrix< DDRMat >            & aResidual,
                                                                      Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolators
            Field_Interpolator* phi = aFieldInterpolators( 0 );
            Field_Interpolator* vN  = aFieldInterpolators( 1 );

            //compute the residual
            aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * phi->gradx( 1 ) );

            // set the jacobian size
            aJacobians.resize( 1 );

            // compute the jacobian Jphiphi
            aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val() * phi->Bx() );

            // compute the jacobian JphivN
            uint tvNNumOfDofs = vN->get_number_of_fields()*vN->get_number_of_space_time_bases();
            aJacobians( 1 ) = trans( phi->N() ) * reshape( trans( vN->N() ) * trans( phi->gradx( 1 ) ), 1, tvNNumOfDofs );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
