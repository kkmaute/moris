
#include "cl_FEM_IWG_Olsson_CLS_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

//        IWG_Olsson_CLS_Bulk::IWG_Olsson_CLS_Bulk( const real aFieldUpperBound,
//                                                  const real aFieldLowerBound,
//                                                  const real aMethodParameter )
        IWG_Olsson_CLS_Bulk::IWG_Olsson_CLS_Bulk()
        {
            //FIXME set field upper and lower bound
            mPhiUB = 1.0;
            mPhiLB = 0.0;

            //FIXME set Olsson CLS epsilon parameter
            mEpsilon = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS2 };

            // set the active dof type
            mMasterDofTypes = {{ MSI::Dof_Type::LS2 },
                               { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY } };
        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_residual( Matrix< DDRMat >                   & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set field interpolators
            Field_Interpolator* phi  = mMasterFI( 0 );
            Field_Interpolator* nPhi = mMasterFI( 1 );

            // compute residual
            aResidual = trans( phi->N() ) * phi->gradt( 1 )
                      - trans( phi->Bx() ) * ( ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) )
                      - mEpsilon  * dot( phi->gradx( 1 ), nPhi->val() ) ) * trans( nPhi->val() );
        }


//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set field interpolators
            Field_Interpolator* phi  = mMasterFI( 0 );
            Field_Interpolator* nPhi = mMasterFI( 1 );

            // set the jacobian size
            aJacobians.resize( 2 );

            //compute the jacobians
            aJacobians( 0 ) = trans( phi->N() )  * phi->Bt()
                            - trans( phi->Bx() ) * ( ( mPhiUB + mPhiLB - 2 * phi->val()( 0 ) ) * trans( nPhi->val() ) * phi->N()
                                                    - mEpsilon * ( trans( nPhi->val() ) * nPhi->val() * phi->Bx() ) );

           // build the global shape functions matrix for vectorial field nPhi
           uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
           uint tNFieldsNPhi = nPhi->get_number_of_fields();
           Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
           for( uint i = 0; i < tNFieldsNPhi; i++ )
           {
               tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
           }


            aJacobians( 1 ) = - trans( phi->Bx() ) *(
                             ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) ) * tNNPhi
                             - mEpsilon * trans( nPhi->val() ) * trans( phi->gradx( 1 ) ) * tNNPhi
                             - mEpsilon * dot( phi->gradx( 1 ), nPhi->val() ) * tNNPhi );

        }

//------------------------------------------------------------------------------

        void IWG_Olsson_CLS_Bulk::compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                                 Matrix< DDRMat >                   & aResidual )
        {
            // check master field interpolators
            this->check_field_interpolators( mtk::Master_Slave::MASTER );

            // set field interpolators
            Field_Interpolator* phi  = mMasterFI( 0 );
            Field_Interpolator* nPhi = mMasterFI( 1 );

            // compute the residual
            aResidual = trans( phi->N() ) * phi->gradt( 1 )
                      - trans( phi->Bx() ) * ( ( phi->val()( 0 ) - mPhiLB ) * (mPhiUB - phi->val()( 0 ) )
                      - mEpsilon  * dot( phi->gradx( 1 ), nPhi->val() ) ) * nPhi->val();

            // set the jacobian size
            aJacobians.resize( 2 );

            //compute the jacobians
            aJacobians( 0 ) = trans( phi->N() )  * phi->Bt()
                             - trans( phi->Bx() ) * ( ( mPhiUB + mPhiLB - 2 * phi->val()( 0 ) ) * trans( nPhi->val() ) * phi->N()
                                                     - mEpsilon * ( trans( nPhi->val() ) * nPhi->val() * phi->Bx() ) );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            aJacobians( 1 ) = - trans( phi->Bx() ) *(
                              ( phi->val()( 0 ) - mPhiLB ) * ( mPhiUB - phi->val()( 0 ) ) * tNNPhi
                              - mEpsilon * trans( nPhi->val() ) * trans( phi->gradx( 1 ) ) * tNNPhi
                              - mEpsilon * dot( phi->gradx( 1 ), nPhi->val() ) * tNNPhi );
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
