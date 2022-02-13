
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        // internal threshold for zero
        const real mEpsilon = 1e-16;

        // internal threshold for wall distance
        const real mWallDistanceEpsilon = 1e-12;

        // Spalart Allmaras model constants
        real mCb1 = 0.1355;
        real mCb2 = 0.6220;
        real mSigma = 2.0/3.0;
        real mKappa = 0.41;
        real mCw1 = mCb1 / std::pow( mKappa, 2.0 ) + ( 1.0 + mCb2 ) / mSigma;
        real mCw2 = 0.3;
        real mCw3 = 2.0;
        real mCt3 = 1.2;
        real mCt4 = 0.5;
        real mCv1 = 7.1;
        real mCv2 = 0.7;
        real mCv3 = 0.9;
        real mRLim = 10.0;
        real mCn1 = 16.0;

        //------------------------------------------------------------------------------

        real compute_chi(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // get the viscosity dof FI
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // compute chi
            return tFIViscosity->val()( 0 ) / aPropKinViscosity->val()( 0 );
        }

        //------------------------------------------------------------------------------

        void compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidu
            adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

            // if derivative dof type is viscosity dof type
            if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
            {
                adchidu += tDerFI->N() / aPropKinViscosity->val()( 0 );
            }

            // if viscosity property depends on derivative dof type
            if( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute chi
                real tChi = compute_chi(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // add contribution to derivative
                adchidu -= tChi * aPropKinViscosity->dPropdDOF( aDofTypes ) /
                        aPropKinViscosity->val()( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void compute_dchidx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & adchidx )
        {
            // get the viscosity dof FI
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dchidx
            adchidx = tFIViscosity->gradx( 1 ) / aPropKinViscosity->val()( 0 );

            // if kinematic viscosity depends on space
            if( aPropKinViscosity->check_space_dependency( 1 ) )
            {
                // add contribution of kinematic viscosity space derivative
                adchidx -= tChi * aPropKinViscosity->dnPropdxn( 1 ) / aPropKinViscosity->val()( 0 );
            }
        }

        //------------------------------------------------------------------------------

        void compute_dchidxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidxdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adchidxdu
            adchidxdu.set_size(
                    tDerFI->gradx( 1 ).n_rows(),
                    tDerFI->get_number_of_space_time_coefficients() );

            // if kinematic viscosity depends on space
            if( aPropKinViscosity->check_space_dependency( 1 ) )
            {
                // compute dchidu
                Matrix< DDRMat > tdChidu;
                compute_dchidu(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdChidu );

                // add contribution of kinematic viscosity space derivative
                adchidxdu = - 1.0 * aPropKinViscosity->dnPropdxn( 1 ) * tdChidu / aPropKinViscosity->val()( 0 );
            }
            else
            {
                adchidxdu.fill( 0.0 );
            }

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // if dof type is viscosity
            if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
            {
                adchidxdu += tDerFI->dnNdxn( 1 ) / aPropKinViscosity->val()( 0 );
            }

            // if viscosity property depends on dof type
            if( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
            {
                // compute chi
                real tChi = compute_chi(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                adchidxdu -= tFIViscosity->gradx( 1 ) *
                        aPropKinViscosity->dPropdDOF( aDofTypes ) / std::pow( aPropKinViscosity->val()( 0 ), 2 );

                // if kinematic viscosity depends on space
                if( aPropKinViscosity->check_space_dependency( 1 ) )
                {
                    adchidxdu += tChi * aPropKinViscosity->dnPropdxn( 1 ) *
                            aPropKinViscosity->dPropdDOF( aDofTypes ) / std::pow( aPropKinViscosity->val()( 0 ), 2 );
                    // FIXME dPropddxdu
                    // - tChi * aPropKinViscosity->dPropdxdu() / aPropKinViscosity->val()
                }
            }
        }

        //------------------------------------------------------------------------------

        real compute_fv1(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // compute chi, chi³
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi³
            real tChi3 = std::pow( tChi, 3.0 );

            // compute fv1
            return tChi3 / ( tChi3 + std::pow( mCv1, 3.0 ) );
        }

        //------------------------------------------------------------------------------

        void compute_dfv1du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1du )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            compute_dchidu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidu );

            // compute adfv1du
            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu /
                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

        //------------------------------------------------------------------------------

        void compute_dfv1dx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & adfv1dx )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            compute_dchidx(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    tdchidx );

            // compute dfv1dx
            adfv1dx =
                    3.0 * mCv1 * std::pow( tChi, 2.0 ) * tdchidx /
                    std::pow( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), 2.0 );
        }

        //------------------------------------------------------------------------------

        void compute_dfv1dxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1dxdu )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi^2
            real tChi2 = std::pow( tChi, 2.0 );

            // compute chi³
            real tChi3 = std::pow( tChi, 3.0 );

            // compute cv1³
            real tCv13 = std::pow( mCv1, 3.0 );

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            compute_dchidx(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    tdchidx );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            compute_dchidu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidu );

            // compute dchidxdu
            Matrix< DDRMat > tdchidxdu;
            compute_dchidxdu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidxdu );

            // compute dfv1dxdu
            adfv1dxdu = 3.0 * mCv1 * (
                    2.0 * tChi * ( tCv13 - 2.0 * tChi3 ) * tdchidx * tdchidu +
                    tChi2 * ( tChi3 + tCv13 ) * tdchidxdu ) /
                            std::pow( tChi3 + tCv13, 3.0 );
        }

        //------------------------------------------------------------------------------

        real compute_fv2(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute fv1
            real tFv1 = compute_fv1(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute fv2
            return 1.0 - tChi / ( 1 + tChi * tFv1 );
        }

        //------------------------------------------------------------------------------

        void compute_dfv2du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv2du )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            compute_dchidu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidu );

            // compute fv1
            real tFv1 = compute_fv1(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dfv1du
            Matrix< DDRMat > tdfv1du;
            compute_dfv1du(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdfv1du );

            // compute adfv2du
            adfv2du = ( std::pow( tChi, 2.0 ) * tdfv1du - tdchidu ) /
                    ( std::pow( 1.0 + tChi * tFv1, 2.0 ) );
        }

    } /* namespace fem */
} /* namespace moris */

