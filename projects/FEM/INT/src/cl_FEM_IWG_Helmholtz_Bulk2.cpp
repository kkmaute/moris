
#include "cl_FEM_IWG_Helmholtz_Bulk2.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"

namespace moris
{
    namespace fem
    {
//------------------------------------------------------------------------------

        IWG_Helmholtz_Bulk2::IWG_Helmholtz_Bulk2()
        {
            //FIXME set the Helmholtz filter parameter
            mFilterParam = 1.0;

            //FIXME set the sharpness parameter for smoothed Heaviside function
            mSharpParam = 1.0;

            //FIXME set the element size
            mHe = 1.0;

            // set the residual dof type
            mResidualDofType = { MSI::Dof_Type::VX };

            // set the active dof type
            mActiveDofTypes = { { MSI::Dof_Type::VX } ,
                                { MSI::Dof_Type::LS1 } };
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_residual( Matrix< DDRMat >            & aResidual,
                                                    Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* vN  = aFieldInterpolators( 0 );
            Field_Interpolator* phi = aFieldInterpolators( 1 );

            //FIXME set unfiltered velocity value
            uint tVNBases = vN->get_number_of_space_time_bases();
            uint tVNFields = vN->get_number_of_fields();
            Matrix< DDRMat > aVHat( tVNBases, tVNFields, 1.0 );

            // compute norm( phi )
            real tNormPhi = norm( phi->gradx( 1 ) );
            if( tNormPhi < 1.0e-12 )
            {
                tNormPhi = 1.0e-12;
            }

            // compute Dirac filter
            real tDiracFilter = 0.0;
            if ( phi->val()( 0 ) < 3 * mHe )
            {
                real tTanh   = std::tanh( mSharpParam * phi->val()( 0 ) / tNormPhi );
                tDiracFilter = 0.5 * mSharpParam * ( 1.0 - std::pow( tTanh, 2 ) );
            }

            // compute the residual
            aResidual = mFilterParam * trans( vN->Bx() ) * vN->gradx( 1 )
                      + trans( vN->N() ) * ( vN->val() - vN->N() * aVHat ) * tDiracFilter;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_jacobian( Cell< Matrix< DDRMat > >    & aJacobians,
                                                    Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* vN  = aFieldInterpolators( 0 );
            Field_Interpolator* phi = aFieldInterpolators( 1 );

            //FIXME set unfiltered velocity value
            uint tVNBases = vN->get_number_of_space_time_bases();
            uint tVNFields = vN->get_number_of_fields();
            Matrix< DDRMat > aVHat( tVNBases, tVNFields, 1.0 );

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

            // compute Dirac filter and derivative wrt phi / norm( phi )
            real tDiracFilter = 0.0;
            real tDDiracFilter = 0.0;
            if ( phi->val()( 0 ) < 3 * mHe )
            {
                real tTanh = std::tanh( mSharpParam * phi->val()( 0 ) / tNormPhi );
                tDiracFilter  = 0.5 * mSharpParam * ( 1.0 - std::pow( tTanh, 2 ) );
                tDDiracFilter = - 2 * mSharpParam * tDiracFilter * tTanh;
            }

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute the jacobian j_vN_vN
            aJacobians( 0 ) = mFilterParam * trans( vN->Bx() ) * vN->Bx() + trans( vN->N() ) * vN->N() * tDiracFilter;

            // compute the jacobian j_vN_phi
            aJacobians( 1 ) = trans( vN->N() ) * ( vN->val() - vN->N() * aVHat )
                            * tDDiracFilter * ( phi->N() * tNormPhi  - phi->val()( 0 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;
        }

//------------------------------------------------------------------------------

        void IWG_Helmholtz_Bulk2::compute_jacobian_and_residual( Cell< Matrix< DDRMat > >    & aJacobians,
                                                                 Matrix< DDRMat >            & aResidual,
                                                                 Cell< Field_Interpolator* > & aFieldInterpolators )
        {
            // set field interpolator
            Field_Interpolator* vN  = aFieldInterpolators( 0 );
            Field_Interpolator* phi = aFieldInterpolators( 1 );

            //FIXME set unfiltered velocity value
            uint tVNBases = vN->get_number_of_space_time_bases();
            uint tVNFields = vN->get_number_of_fields();
            Matrix< DDRMat > aVHat( tVNBases, tVNFields, 1.0 );

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

            // compute Dirac filter and derivative wrt phi / norm( phi )
            real tDiracFilter = 0.0;
            real tDDiracFilter = 0.0;
            if ( phi->val()( 0 ) < 3 * mHe )
            {
                real tTanh = std::tanh( mSharpParam * phi->val()( 0 ) / tNormPhi );
                tDiracFilter  = 0.5 * mSharpParam * ( 1.0 - std::pow( tTanh, 2 ) );
                tDDiracFilter = - 2 * mSharpParam * tDiracFilter * tTanh;
            }

            // compute the residual r_vN
            aResidual = mFilterParam * trans( vN->Bx() ) * vN->gradx( 1 )
                      + trans( vN->N() ) * ( vN->val() - vN->N() * aVHat ) * tDiracFilter;

            // set the jacobian size
            aJacobians.resize( 2 );

            // compute the jacobian j_vN_vN
            aJacobians( 0 ) = mFilterParam * trans( vN->Bx() ) * vN->Bx() + trans( vN->N() ) * vN->N() * tDiracFilter;

           // compute the jacobian j_vN_phi
           aJacobians( 1 ) = trans( vN->N() ) * ( vN->val() - vN->N() * aVHat )
                           * tDDiracFilter * ( phi->N() * tNormPhi  - phi->val()( 0 ) * trans( tDNormPhiDPhiHat ) ) / std::pow( tNormPhi, 2 ) ;
        }

//------------------------------------------------------------------------------
    } /* namespace fem */
} /* namespace moris */
