
#include "cl_FEM_IWG_Hamilton_Jacobi_Bulk2.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"

namespace moris
{
    namespace fem
    {

//------------------------------------------------------------------------------

        IWG_Hamilton_Jacobi_Bulk2::IWG_Hamilton_Jacobi_Bulk2()
        {
            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::LS1 };

            // set the active dof type
            mMasterDofTypes = { { MSI::Dof_Type::LS1 },
                                { MSI::Dof_Type::VX } };
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk2::compute_residual( Matrix< DDRMat >                   & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->gradx( 1 ) ) * phi->Bx() / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
            }

           //compute the residual
           aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * dot( phi->gradx( 1 ), phi->gradx( 1 ) ) / tNormPhi );
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk2::compute_jacobian( moris::Cell< Matrix< DDRMat > >    & aJacobians )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

            // compute norm( phi ) and derivative wrt phiHat
            real tNormPhi                     = norm( phi->gradx( 1 ) );
            Matrix< DDRMat > tDNormPhiDPhiHat = trans( phi->Bx() ) * phi->gradx( 1 ) / tNormPhi;

            // If all values of level set in this element are the same,
            // then gradient is zero, protect from going to NAN/inf
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
                uint tNPhiCoeff = phi->get_number_of_space_time_coefficients();
                tDNormPhiDPhiHat.set_size( tNPhiCoeff, 1, 0.0 );
            }

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute the jacobian Jphiphi
            aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val()( 0 ) * trans( tDNormPhiDPhiHat ) );

            // compute the jacobian JphivN
            aJacobians( 1 ) = trans( phi->N() ) * vN->N() * dot( phi->gradx( 1 ), phi->gradx( 1 ) ) / tNormPhi;
        }

//------------------------------------------------------------------------------

        void IWG_Hamilton_Jacobi_Bulk2::compute_jacobian_and_residual( moris::Cell< Matrix< DDRMat > >    & aJacobians,
                                                                       Matrix< DDRMat >                   & aResidual )
        {
            // set field interpolators
            Field_Interpolator* phi = mMasterFI( 0 );
            Field_Interpolator* vN  = mMasterFI( 1 );

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

            //compute the residual
            aResidual = trans( phi->N() ) * ( phi->gradt( 1 ) + vN->val() * dot( phi->gradx( 1 ), phi->gradx( 1 ) ) / tNormPhi );

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute the jacobian Jphiphi
            aJacobians( 0 ) = trans( phi->N() ) * ( phi->Bt() + vN->val()( 0 ) * trans( tDNormPhiDPhiHat ) );

            // compute the jacobian JphivN
            aJacobians( 1 ) = trans( phi->N() ) * vN->N() * dot( phi->gradx( 1 ), phi->gradx( 1 ) ) / tNormPhi;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
