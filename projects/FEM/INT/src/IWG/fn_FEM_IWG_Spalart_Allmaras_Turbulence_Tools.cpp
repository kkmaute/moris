/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.cpp
 *
 */

#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"
#include "fn_clip_value.hpp"

namespace moris::fem
{
    // Spalart Allmaras model constants
    const real mCv1 = 7.1;

    const real mEpsilon      = 1e-10;    // needs to be consistent with threshold set in
    const real mEpsilonDeriv = 1e-60;    // cl_FEM_CM_Spalart_Allmaras_Turbulence.hpp

    //------------------------------------------------------------------------------

    real
    compute_chi(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity )
    {
        // get the viscosity dof FI
        Field_Interpolator* tFIViscosity =
                aLeaderFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

        // compute chi
        return tFIViscosity->val()( 0 ) / aPropKinViscosity->val()( 0 );
    }

    //------------------------------------------------------------------------------

    void
    compute_dchidu(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            const Vector< MSI::Dof_Type >&     aDofTypes,
            Matrix< DDRMat >&                  adchidu )
    {
        // get the derivative dof FIs
        Field_Interpolator* tDerFI =
                aLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init adchidu
        adchidu.set_size( 1, tDerFI->get_number_of_space_time_coefficients(), 0.0 );

        // if derivative dof type is viscosity dof type
        if ( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
        {
            adchidu += tDerFI->N() / aPropKinViscosity->val()( 0 );
        }

        // if viscosity property depends on derivative dof type
        if ( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aLeaderFIManager,
                    aPropKinViscosity );

            // add contribution to derivative
            adchidu -=
                    tChi * aPropKinViscosity->dPropdDOF( aDofTypes ) / aPropKinViscosity->val()( 0 );
        }
    }

    //------------------------------------------------------------------------------

    void
    compute_dchidx(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            Matrix< DDRMat >&                  adchidx )
    {
        // get the viscosity dof FI
        Field_Interpolator* tFIViscosity =
                aLeaderFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

        // compute dchidx
        adchidx = tFIViscosity->gradx( 1 ) / aPropKinViscosity->val()( 0 );

        // if kinematic viscosity depends on space
        if ( aPropKinViscosity->check_space_dependency() )
        {
            // assume that kinematic viscosity prop does not depend on x
            MORIS_ERROR( false,           //
                    "compute_dchidx -"    //
                    "Dependence of kinematic viscosity on space not accounted for." );
        }
    }

    //------------------------------------------------------------------------------

    void
    compute_dchidxdu(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            const Vector< MSI::Dof_Type >&     aDofTypes,
            Matrix< DDRMat >&                  adchidxdu )
    {
        // get the derivative dof FIs
        Field_Interpolator* tDerFI =
                aLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init adchidxdu
        adchidxdu.set_size(
                tDerFI->gradx( 1 ).n_rows(),
                tDerFI->get_number_of_space_time_coefficients() );

        // if kinematic viscosity depends on space
        if ( aPropKinViscosity->check_space_dependency() )
        {
            // assume that kinematic viscosity prop does not depend on x
            MORIS_ERROR( false,             //
                    "compute_dchidxdu -"    //
                    "Dependence of kinematic viscosity on space not accounted for." );
        }
        else
        {
            adchidxdu.fill( 0.0 );
        }

        // get the residual dof FI (here viscosity)
        Field_Interpolator* tFIViscosity =
                aLeaderFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

        // if dof type is viscosity
        if ( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
        {
            adchidxdu += tDerFI->dnNdxn( 1 ) / aPropKinViscosity->val()( 0 );
        }

        // if viscosity property depends on dof type
        if ( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
        {
            adchidxdu -=
                    tFIViscosity->gradx( 1 ) * aPropKinViscosity->dPropdDOF( aDofTypes )
                    / std::pow( aPropKinViscosity->val()( 0 ), 2 );

            // if kinematic viscosity depends on space
            if ( aPropKinViscosity->check_space_dependency() )
            {
                // assume that kinematic viscosity prop does not depend on x
                MORIS_ERROR( false,              //
                        "compute_dchidxdu - "    //
                        "Dependence of kinematic viscosity on space not accounted for." );
            }
        }
    }

    //------------------------------------------------------------------------------

    real
    compute_fv1(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity )
    {
        // compute chi, chi³
        real tChi = compute_chi(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity );

        // compute chi³
        real tChi3 = std::pow( tChi, 3.0 );

        // threshold deno (for consistency with derivative computation)
        real tFv1Deno = clip_value( tChi3 + std::pow( mCv1, 3.0 ), mEpsilon );

        // compute fv1
        return tChi3 / tFv1Deno;
    }

    //------------------------------------------------------------------------------

    void
    compute_dfv1du(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            const Vector< MSI::Dof_Type >&     aDofTypes,
            Matrix< DDRMat >&                  adfv1du )
    {
        // compute chi
        real tChi = compute_chi(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity );

        // compute dchidu
        Matrix< DDRMat > tdchidu;
        compute_dchidu(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity,
                aDofTypes,
                tdchidu );

        // threshold deno (for consistency with derivative computation)
        real tFv1Deno = clip_value( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), mEpsilon );

        if ( std::abs( tFv1Deno ) > mEpsilon )
        {
            // threshold denominator
            real tFv1Deno2 = clip_value( std::pow( tFv1Deno, 2.0 ), mEpsilonDeriv );

            // compute adfv1du
            adfv1du = 3.0 * std::pow( mCv1, 3.0 ) * std::pow( tChi, 2.0 ) * tdchidu / tFv1Deno2;
        }
        else
        {
            adfv1du = 3.0 * std::pow( tChi, 2.0 ) * tdchidu / tFv1Deno;
        }
    }

    //------------------------------------------------------------------------------

    void
    compute_dfv1dx(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            Matrix< DDRMat >&                  adfv1dx )
    {
        // compute chi
        real tChi = compute_chi(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity );

        // compute dchidx
        Matrix< DDRMat > tdchidx;
        compute_dchidx(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity,
                tdchidx );

        // threshold deno (for consistency with derivative computation)
        real tFv1Deno = clip_value( std::pow( tChi, 3.0 ) + std::pow( mCv1, 3.0 ), mEpsilon );

        if ( std::abs( tFv1Deno ) > mEpsilon )
        {
            // threshold denominator
            real tFv1Deno2 = clip_value( std::pow( tFv1Deno, 2.0 ), mEpsilonDeriv );

            // compute dfv1dx
            adfv1dx = 3.0 * mCv1 * std::pow( tChi, 2.0 ) * tdchidx / tFv1Deno2;
        }
        else
        {
            // compute dfv1dx
            adfv1dx = 3.0 * std::pow( tChi, 2.0 ) * tdchidx / tFv1Deno;
        }
    }

    //------------------------------------------------------------------------------

    void
    compute_dfv1dxdu(
            const Vector< MSI::Dof_Type >&     aViscosityDofGroup,
            Field_Interpolator_Manager*        aLeaderFIManager,
            const std::shared_ptr< Property >& aPropKinViscosity,
            const Vector< MSI::Dof_Type >&     aDofTypes,
            Matrix< DDRMat >&                  adfv1dxdu )
    {
        // compute chi
        real tChi = compute_chi(
                aViscosityDofGroup,
                aLeaderFIManager,
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
                aLeaderFIManager,
                aPropKinViscosity,
                tdchidx );

        // compute dchidu
        Matrix< DDRMat > tdchidu;
        compute_dchidu(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity,
                aDofTypes,
                tdchidu );

        // compute dchidxdu
        Matrix< DDRMat > tdchidxdu;
        compute_dchidxdu(
                aViscosityDofGroup,
                aLeaderFIManager,
                aPropKinViscosity,
                aDofTypes,
                tdchidxdu );

        // threshold deno (for consistency with derivative computation)
        real tFv1Deno = clip_value( tChi3 + tCv13, mEpsilon );

        if ( std::abs( tFv1Deno ) > mEpsilon )
        {
            // threshold denominator
            real tFv1Deno3 = clip_value( std::pow( tFv1Deno, 3.0 ), mEpsilonDeriv );

            // compute dfv1dxdu
            adfv1dxdu =
                    3.0 * mCv1 * ( 2.0 * tChi * ( tCv13 - 2.0 * tChi3 ) * tdchidx * tdchidu + tChi2 * ( tChi3 + tCv13 ) * tdchidxdu ) / tFv1Deno3;
        }
        else
        {
            // compute dfv1dxdu
            adfv1dxdu = 3.0 * tChi * ( 2.0 * tdchidx * tdchidu + tChi * tdchidxdu ) / tFv1Deno;
        }
    }

}    // namespace moris::fem
