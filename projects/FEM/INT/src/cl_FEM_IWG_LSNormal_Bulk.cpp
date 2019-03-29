
#include "cl_FEM_IWG_LSNormal_Bulk.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "op_div.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        IWG_LSNormal_Bulk::IWG_LSNormal_Bulk()
        {
            // set the residual dof type
            //FIXME should be NLSX, NLSY, NLSZ
            //mResidualDofType = { MSI::Dof_Type::NLSX };
            mResidualDofType = { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY };
            //mResidualDofType = { MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ };

            // set the active dof type
            //mActiveDofTypes = {{ MSI::Dof_Type::NLSX },
            //                   { MSI::Dof_Type::LS1 }};
            mActiveDofTypes = {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY },
                               { MSI::Dof_Type::LS1 }};
            //mActiveDofTypes = {{ MSI::Dof_Type::NLSX, MSI::Dof_Type::NLSY, MSI::Dof_Type::NLSZ },
            //                   { MSI::Dof_Type::LS1 }};
        }

//------------------------------------------------------------------------------

        void IWG_LSNormal_Bulk::compute_residual( Matrix< DDRMat >            & aResidual,
                                                  Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolators
            Field_Interpolator* nPhi = aFieldInterpolators( 0 );
            Field_Interpolator* phi  = aFieldInterpolators( 1 );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            // compute norm( phi )
            real tNormPhi = norm( phi->gradx( 1 ) );
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
            }

            // compute residual
            aResidual = trans( tNNPhi ) * ( trans( nPhi->val() ) - phi->gradx( 1 ) / tNormPhi );
        }

//------------------------------------------------------------------------------

        void IWG_LSNormal_Bulk::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                  Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolators
            Field_Interpolator* nPhi = aFieldInterpolators( 0 );
            Field_Interpolator* phi  = aFieldInterpolators( 1 );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->Bx() ) * phi->gradx( 1 ) / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
                uint tNPhiBases = phi->get_number_of_space_time_bases();
                tDNormPhiDPhiHat.set_size( tNPhiBases, 1, 0.0 );
            }

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute j_nPhi_nPhi
            aJacobians( 0 ) = trans( tNNPhi ) * tNNPhi;

            // compute j_nPhi_phi
            aJacobians( 1 ) = - trans( tNNPhi ) * ( phi->Bx() * tNormPhi - phi->gradx( 1 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;
        }

//------------------------------------------------------------------------------

        void IWG_LSNormal_Bulk::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                               Matrix< DDRMat >            & aResidual,
                                                               Cell< Field_Interpolator* > & aFieldInterpolators)
        {
            // set field interpolators
            Field_Interpolator* nPhi = aFieldInterpolators( 0 );
            Field_Interpolator* phi  = aFieldInterpolators( 1 );

            // build the global shape functions matrix for vectorial field nPhi
            uint tNBasesNPhi  = nPhi->get_number_of_space_time_bases();
            uint tNFieldsNPhi = nPhi->get_number_of_fields();
            Matrix< DDRMat > tNNPhi( tNFieldsNPhi, tNFieldsNPhi * tNBasesNPhi, 0.0 );
            for( uint i = 0; i < tNFieldsNPhi; i++ )
            {
                tNNPhi({i,i},{i * tNBasesNPhi, (i+1) * tNBasesNPhi - 1}) = nPhi->N().get_row( 0 );
            }

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->gradx( 1 ) ) * phi->Bx() / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
                uint tNPhiBases = phi->get_number_of_space_time_bases();
                tDNormPhiDPhiHat.set_size( tNPhiBases, 1, 0.0 );
            }

            // compute residual
            aResidual = trans( tNNPhi ) * ( trans( nPhi->val() ) - phi->gradx( 1 ) / tNormPhi );

            // compute j_nPhi_nPhi
            aJacobians( 0 ) = trans( tNNPhi ) * tNNPhi;

            // compute j_nPhi_phi
            aJacobians( 1 ) = - trans( tNNPhi ) * ( phi->Bx() * tNormPhi - phi->gradx( 1 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;

        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
