
#include "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        // internal threshold for zero
        const real mEpsilon = 1e-18;

        // internal threshold for wall distance
        const real mWallDistanceEpsilon = 1e-6;

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

        real compute_production_term(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // init production term
            real tProductionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity >= 0.0 )
            {
                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute Stilde
                real tSTilde = compute_stilde(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute production term
                tProductionTerm = mCb1 * ( 1.0 - tFt2 ) * tSTilde * tModViscosity;
            }
            // if viscosity is negative
            else
            {
                // compute s
                real tS = compute_s(
                        aVelocityDofGroup,
                        aMasterFIManager );

                // compute production term
                tProductionTerm = mCb1 * ( 1.0 - mCt3 ) * tS * tModViscosity;
            }

            return tProductionTerm;
        }

        //------------------------------------------------------------------------------

        void compute_dproductiontermdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adproductiondu )
        {
            //
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adproductiondu
            adproductiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dft2du
                Matrix< DDRMat > tdft2du;
                compute_dft2du(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdft2du );

                // compute STilde
                real tSTilde = compute_stilde(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute dSTildedu
                Matrix< DDRMat > tdstildedu;
                compute_dstildedu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdstildedu );

                // add contribution to dproductiondu
                adproductiondu =
                        - mCb1 * tSTilde * tModViscosity * tdft2du +
                        mCb1 * ( 1 - tFt2 ) * tModViscosity * tdstildedu;

                // if derivative dof type is viscosity dof type
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dproductiondu
                    adproductiondu +=
                            mCb1 * ( 1.0 - tFt2 ) * tSTilde * tFIModViscosity->N();
                }
            }
            // if viscosity is negative
            else
            {
                // compute s
                real tS = compute_s(
                        aVelocityDofGroup,
                        aMasterFIManager );

                // compute dsdu
                Matrix< DDRMat > tdsdu;
                compute_dsdu(
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aDofTypes,
                        tdsdu );

                // add contribution to dproductiondu
                adproductiondu =
                        mCb1 * ( 1.0 - mCt3 ) * tModViscosity * tdsdu;

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dproductiondu
                    adproductiondu +=
                            mCb1 * ( 1.0 - mCt3 ) * tS * tFIModViscosity->N();
                }
            }
        }

        //------------------------------------------------------------------------------

        real compute_production_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // init production coefficient
            real tProductionCoeff = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute Stilde
                real tSTilde = compute_stilde(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute production coefficient
                tProductionCoeff = mCb1 * ( 1.0 - tFt2 ) * tSTilde;
            }
            // if viscosity is negative
            else
            {
                // compute s
                real tS = compute_s(
                        aVelocityDofGroup,
                        aMasterFIManager );

                // compute production coefficient
                tProductionCoeff = mCb1 * ( 1.0 - mCt3 ) * tS;
            }

            return tProductionCoeff;
        }

        //------------------------------------------------------------------------------

        void compute_dproductioncoefficientdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adproductiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adproductiondu
            adproductiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dft2du
                Matrix< DDRMat > tdft2du;
                compute_dft2du(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity ,
                        aDofTypes,
                        tdft2du );

                // compute Stilde
                real tSTilde = compute_stilde(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute dSTildedu
                Matrix< DDRMat > tdstildedu;
                compute_dstildedu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdstildedu );

                // add contribution to dproductiondu
                adproductiondu = - mCb1 * tSTilde * tdft2du +
                        mCb1 * ( 1 - tFt2 ) * tdstildedu;
            }
            // if viscosity is negative
            else
            {
                // compute dsdu
                Matrix< DDRMat > tdsdu;
                compute_dsdu(
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aDofTypes,
                        tdsdu );

                // add contribution to dproductiondu
                adproductiondu = mCb1 * ( 1.0 - mCt3 ) * tdsdu;
            }
        }

        //------------------------------------------------------------------------------

        real compute_wall_destruction_term(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // init wall destruction term
            real tWallDestructionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute fw
                real tFw = compute_fw(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute wall destruction term
                tWallDestructionTerm =
                        ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                        std::pow( tModViscosity / tWallDistance, 2.0 );
            }
            // if viscosity is negative
            else
            {
                // compute wall destruction term
                tWallDestructionTerm = - mCw1 * std::pow( tModViscosity / tWallDistance, 2.0 );
            }

            return tWallDestructionTerm;
        }

        //------------------------------------------------------------------------------

        void compute_dwalldestructiontermdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwalldestructiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adwalldestructiondu
            adwalldestructiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients());

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute fw
                real tFw = compute_fw(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute dfwfu
                Matrix< DDRMat > tdfwdu;
                compute_dfwdu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdfwdu );

                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dft2fu
                Matrix< DDRMat > tdft2du;
                compute_dft2du(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdft2du );

                // add contribution from ft2 and fw to dwalldestructiondu
                adwalldestructiondu =
                        ( mCw1 * tdfwdu - mCb1 * tdft2du / std::pow( mKappa, 2.0 ) ) *
                        std::pow( tModViscosity / tWallDistance, 2.0 );

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu +=
                            ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                            2.0 * tModViscosity * tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }

                // if wall distance depends on derivative dof type
                if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                        ( tWallDistance > mWallDistanceEpsilon ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu -=
                            2.0 * ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                            std::pow( tModViscosity, 2.0 ) *
                            aPropWallDistance->dPropdDOF( aDofTypes ) /
                            std::pow( tWallDistance, 3.0 );
                }
            }
            // if viscosity is negative
            else
            {
                // fill with zeros
                adwalldestructiondu.fill( 0.0 );

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu -=
                            2.0 * mCw1 * tModViscosity * tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }

                // if wall distance depends on derivative dof type
                if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                        ( tWallDistance > mWallDistanceEpsilon ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu +=
                            2.0 * mCw1 * std::pow( tModViscosity, 2.0 ) *
                            aPropWallDistance->dPropdDOF( aDofTypes ) /
                            std::pow( tWallDistance, 3.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        real compute_wall_destruction_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // init wall destruction coefficient
            real tWallDestructionCoeff = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute fw
                real tFw = compute_fw(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute wall destruction coefficient
                tWallDestructionCoeff =
                        ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                        tModViscosity / std::pow( tWallDistance, 2.0 );
            }
            // if viscosity is negative
            else
            {
                // compute wall destruction coefficient
                tWallDestructionCoeff = - mCw1 * tModViscosity / std::pow( tWallDistance, 2.0 );
            }

            return tWallDestructionCoeff;
        }

        //------------------------------------------------------------------------------

        void compute_dwalldestructioncoefficientdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwalldestructiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init adwalldestructiondu
            adwalldestructiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute fw
                real tFw = compute_fw(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // compute dfwfu
                Matrix< DDRMat > tdfwdu;
                compute_dfwdu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdfwdu );

                // compute ft2
                real tFt2 = compute_ft2(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dft2fu
                Matrix< DDRMat > tdft2du;
                compute_dft2du(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdft2du );

                // add contribution to dwalldestructiondu
                adwalldestructiondu =
                        ( mCw1 * tdfwdu - mCb1 * tdft2du / std::pow( mKappa, 2.0 ) ) *
                        tModViscosity / std::pow( tWallDistance, 2.0 );

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu +=
                            ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                            tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }

                // if wall distance depends on derivative dof type
                if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                        ( tWallDistance > mWallDistanceEpsilon ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu -=
                            2.0 * ( mCw1 * tFw - mCb1 * tFt2 / std::pow( mKappa, 2.0 ) ) *
                            std::pow( tModViscosity, 2.0 ) *
                            aPropWallDistance->dPropdDOF( aDofTypes ) /
                            std::pow( tWallDistance, 3.0 );
                }
            }
            // if viscosity is negative
            else
            {
                // fill with zeros
                adwalldestructiondu.fill( 0.0 );

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu -=
                            mCw1 * tFIModViscosity->N() /
                            std::pow( tWallDistance, 2.0 );
                }

                // if wall distance depends on derivative dof type
                if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                        ( tWallDistance > mWallDistanceEpsilon ) )
                {
                    // add contribution to dwalldestructiondu
                    adwalldestructiondu +=
                            2.0 * mCw1 * std::pow( tModViscosity, 2.0 ) *
                            aPropWallDistance->dPropdDOF( aDofTypes ) /
                            std::pow( tWallDistance, 3.0 );
                }
            }
        }

        //------------------------------------------------------------------------------

        real compute_diffusion_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // init diffusion coeff
            real tDiffusionTerm = 0.0;

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // get the fluid kinematic viscosity value
            real tKinViscosity = aPropKinViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity ) / mSigma;
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = compute_fn(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute diffusion term
                tDiffusionTerm = ( tKinViscosity + tModViscosity * tFn ) / mSigma;
            }

            return tDiffusionTerm;
        }

        //------------------------------------------------------------------------------

        void compute_ddiffusiondu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // init ddiffusiondu
            addiffusiondu.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // fill with zeros
                addiffusiondu.fill( 0.0 );

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += aPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = compute_fn(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dfndu
                Matrix< DDRMat > tdfndu;
                compute_dfndu(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdfndu );

                // add contribution from fn to ddiffusiondu
                addiffusiondu = tModViscosity * tdfndu / mSigma;

                // if derivative dof type is viscosity
                if( tDerDofType == aViscosityDofGroup( 0 ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += tFn * tFIModViscosity->N() / mSigma;
                }

                // if kinematic viscosity depends on derivative dof type
                if( aPropKinViscosity->check_dof_dependency( aDofTypes ) )
                {
                    // add contribution to ddiffusiondu
                    addiffusiondu += aPropKinViscosity->dPropdDOF( aDofTypes ) / mSigma;
                }
            }
        }

        //------------------------------------------------------------------------------

        void compute_ddiffusiondx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & addiffusiondx )
        {
            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // init diffusion coeff
            addiffusiondx.set_size( tFIModViscosity->gradx( 1 ).n_rows(), 1 );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // compute diffusion term
                addiffusiondx = tFIModViscosity->gradx( 1 ) / mSigma;

                // if kinematic viscosity depends on space
                if( aPropKinViscosity->check_space_dependency( 1 ) )
                {
                    addiffusiondx += aPropKinViscosity->dnPropdxn( 1 ) / mSigma;
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = compute_fn(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dfndx
                Matrix< DDRMat > dfndx;
                compute_dfndx(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        dfndx );

                // compute diffusion term
                addiffusiondx = (
                        tFIModViscosity->gradx( 1 ) * tFn +
                        tFIModViscosity->val( )( 0 ) * dfndx ) / mSigma;

                // if kinematic viscosity depends on space
                if( aPropKinViscosity->check_space_dependency( 1 ) )
                {
                    addiffusiondx += aPropKinViscosity->dnPropdxn( 1 ) / mSigma;
                }
            }
        }

        //------------------------------------------------------------------------------

        void compute_ddiffusiondxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondxdu )
        {
            // get derivative dof type
            MSI::Dof_Type tDerDofType = aDofTypes( 0 );

            // get the derivative dof FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( tDerDofType );

            // get the viscosity FI
            Field_Interpolator * tFIModViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // init ddiffusiondxdu
            addiffusiondxdu.set_size(
                    tFIModViscosity->gradx( 1 ).n_rows(),
                    tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity value
            real tModViscosity = tFIModViscosity->val()( 0 );

            // if viscosity is positive or zero
            if( tModViscosity > 0.0 )
            {
                // fill with zeros
                addiffusiondxdu.fill( 0.0 );

                // if derivative dof type is viscosity
                if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
                {
                    // compute diffusion term
                    addiffusiondxdu += ( tFIModViscosity->dnNdxn( 1 ) ) / mSigma;

//                    //FIXME missing term aPropKinViscosity->dPropdxdu()
//                    // if kinematic viscosity depends on space
//                    if( aPropKinViscosity->check_space_dependency( 1 ) )
//                    {
//                        addiffusiondxdu += aPropKinViscosity->dPropdxdu / mSigma;
//                    }
                }
            }
            // if viscosity is negative
            else
            {
                // compute fn
                real tFn = compute_fn(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity );

                // compute dfndu
                Matrix< DDRMat > tdfndu;
                compute_dfndu(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdfndu );

                // compute dfndx
                Matrix< DDRMat > dfndx;
                compute_dfndx(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        dfndx );

                // compute dfndx
                Matrix< DDRMat > tdfndxdu;
                compute_dfndxdu(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aDofTypes,
                        tdfndxdu );

                // compute diffusion term
                addiffusiondxdu = ( tFIModViscosity->gradx( 1 ) * tdfndu +
                        tFIModViscosity->val()( 0 ) * tdfndxdu ) / mSigma;

                //FIXME missing term aPropKinViscosity->dPropdxdu()
//                    // if kinematic viscosity depends on space
//                    if( aPropKinViscosity->check_space_dependency( 1 ) )
//                    {
//                        addiffusiondxdu += aPropKinViscosity->dPropdxdu() / mSigma;
//                    }

                // if derivative dof type is viscosity
                if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
                {
                    addiffusiondxdu += ( tFIModViscosity->dnNdxn( 1 ) * tFn +
                            dfndx * tFIModViscosity->N() ) / mSigma;
                }
            }
        }

        //------------------------------------------------------------------------------

        void compute_wij(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                Matrix< DDRMat >                   & aWij )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    aMasterFIManager->get_field_interpolators_for_type( aVelocityDofGroup( 0 ) );

            // get gradient of velocity
            const Matrix< DDRMat > & tGradVelocity = tFIVelocity->gradx( 1 );

            // switch on space dim
            uint tSpaceDim = tFIVelocity->val().numel();
            switch ( tSpaceDim )
            {
                case 2 :
                {
                    // init aWij = [ w11 w12 w21 w22]
                    aWij.set_size( 4, 1, 0.0 );

                    // compute Wij
                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    aWij( 2 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    break;
                }
                case 3 :
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    aWij.set_size( 9, 1, 0.0 );

                    // compute Wij
                    aWij( 1 ) = 0.5 * ( tGradVelocity( 1, 0 ) - tGradVelocity( 0, 1 ) );
                    aWij( 2 ) = 0.5 * ( tGradVelocity( 2, 0 ) - tGradVelocity( 0, 2 ) );
                    aWij( 3 ) = 0.5 * ( tGradVelocity( 0, 1 ) - tGradVelocity( 1, 0 ) );
                    aWij( 5 ) = 0.5 * ( tGradVelocity( 2, 1 ) - tGradVelocity( 1, 2 ) );
                    aWij( 6 ) = 0.5 * ( tGradVelocity( 0, 2 ) - tGradVelocity( 2, 0 ) );
                    aWij( 7 ) = 0.5 * ( tGradVelocity( 1, 2 ) - tGradVelocity( 2, 1 ) );
                    break;
                }
                default:
                    MORIS_ERROR( false, "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools::compute_wij - space dim can only be 2 or 3" );
            }
        }

        //------------------------------------------------------------------------------

        void compute_dwijdu(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwijdu )
        {
            // get the der FI
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =
                    aMasterFIManager->get_field_interpolators_for_type( aVelocityDofGroup( 0 ) );

            // switch on space dim
            uint tSpaceDim = tFIVelocity->val().numel();
            switch ( tSpaceDim )
            {
                case 2 :
                {
                    // init aWij = [ w11 w12 w21 w22]
                    adwijdu.set_size( 4, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if( aDofTypes( 0 ) == aVelocityDofGroup( 0 ) )
                    {
                        // get gradient of velocity
                        const Matrix< DDRMat > & tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )             =   0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );
                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )             = - 0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );
                    }
                    break;
                }
                case 3 :
                {
                    // init aWij = [ w11 w12 w13 w21 w22 w23 w31 w32 w33 ]
                    adwijdu.set_size( 9, tFIDer->get_number_of_space_time_coefficients(), 0.0 );

                    if( aDofTypes( 0 ) == aVelocityDofGroup( 0 ) )
                    {
                        // get gradient of velocity
                        const Matrix< DDRMat > & tdNdxVelocity = tFIVelocity->dnNdxn( 1 );

                        // get number of bases for displacement
                        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

                        // compute adwijdu
                        adwijdu( { 1, 1 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 2, 2 }, { 0, tNumBases - 1 } )                 =   0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 3, 3 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 1 );
                        adwijdu( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 5, 5 }, { tNumBases, 2 * tNumBases - 1 } )     =   0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 5, 5 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) = - 0.5 * tdNdxVelocity.get_row( 1 );

                        adwijdu( { 6, 6 }, { 0, tNumBases - 1 } )                 = - 0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 6, 6 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 0 );

                        adwijdu( { 7, 7 }, { tNumBases, 2 * tNumBases - 1 } )     = - 0.5 * tdNdxVelocity.get_row( 2 );
                        adwijdu( { 7, 7 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =   0.5 * tdNdxVelocity.get_row( 1 );
                    }
                    break;
                }
                default:
                    MORIS_ERROR( false, "fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools::compute_dwijdu - space dim can only be 2 or 3" );
            }
        }

        //------------------------------------------------------------------------------

        real compute_s(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager )
        {
            // compute wij
            Matrix< DDRMat > tWij;
            compute_wij(
                    aVelocityDofGroup,
                    aMasterFIManager,
                    tWij );

            // compute WijWij
            Matrix< DDRMat > tWijWij = trans( tWij ) * tWij;

            // compute s
            real tS = std::sqrt( 2.0 * tWijWij( 0 ) );

            // threshold s for consistency with derivative
            return std::max( tS, mEpsilon );
        }

        //------------------------------------------------------------------------------

        void compute_dsdu(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager ,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsdu
            adsdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients() );

            // compute sbar
            real tS = compute_s(
                    aVelocityDofGroup,
                    aMasterFIManager );

            // if s is greater than threshold
            if( tS > mEpsilon )
            {
                // compute wij
                Matrix< DDRMat > tWij;
                compute_wij(
                        aVelocityDofGroup,
                        aMasterFIManager,
                        tWij );

                // compute dwijdu
                Matrix< DDRMat > tdWijdu;
                compute_dwijdu(
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aDofTypes,
                        tdWijdu );

                // compute dsdu
                adsdu = 2.0 * trans( tWij ) * tdWijdu / tS;
            }
            else
            {
                adsdu.fill( 0.0 );
            }
        }

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

        //------------------------------------------------------------------------------

        real compute_fn(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi³
            real tChi3 = std::pow( tChi, 3 );

            // compute fn
            return ( mCn1 + tChi3 ) / ( mCn1 - tChi3 );
        }

        //------------------------------------------------------------------------------

        void compute_dfndu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adfndxdu
            adfndu.set_size(
                    1,
                    tDerFI->get_number_of_space_time_coefficients() );

            // compute chi, chi², chi³
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi²
            real tChi2 = std::pow( tChi, 2 );

            // compute chi³
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            compute_dchidu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidu );

            // compute adfndu
            adfndu = 6.0 * mCn1 * tChi2 * tdchidu / std::pow( mCn1 - tChi3, 2.0 );
        }

        //------------------------------------------------------------------------------

        void compute_dfndx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & adfndx )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi²
            real tChi2 = std::pow( tChi, 2 );

            // compute chi³
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            compute_dchidx(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    tdchidx );

            // compute adfndx
            adfndx = 6.0 * mCn1 * tChi2 * tdchidx / std::pow( mCn1 - tChi3, 2.0 );
        }

        //------------------------------------------------------------------------------

        void compute_dfndxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndxdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adfndxdu
            adfndxdu.set_size(
                    tDerFI->gradx( 1 ).n_rows(),
                    tDerFI->get_number_of_space_time_coefficients() );

            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute chi²
            real tChi2 = std::pow( tChi, 2 );

            // compute chi³
            real tChi3 = std::pow( tChi, 3 );

            // compute dchidu
            Matrix< DDRMat > tdchidu;
            compute_dchidu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidu );

            // compute dchidx
            Matrix< DDRMat > tdchidx;
            compute_dchidx(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    tdchidx );

            // compute dchidxdu
            Matrix< DDRMat > tdchidxdu;
            compute_dchidxdu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdchidxdu );

            // compute adfndx
            adfndxdu = 6.0 * mCn1 * ( 2.0 * tChi * ( mCn1 + 2.0 * tChi3 ) * tdchidx * tdchidu +
                    ( mCn1 - tChi3 ) * tChi2 * tdchidxdu ) / std::pow( mCn1 - tChi3, 3.0 );
        }

        //------------------------------------------------------------------------------

        real compute_sbar(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // get the viscosity FI
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get wall distance
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // compute fv2
            real tFv2 = compute_fv2(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute s
            return tFv2 * tFIViscosity->val()( 0 ) / std::pow( mKappa * tWallDistance, 2.0 );
        }

        //------------------------------------------------------------------------------

        void compute_dsbardu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsbardu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tFIDer =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init dsbardu
            adsbardu.set_size( 1, tFIDer->get_number_of_space_time_coefficients() );

            // get the viscosity FI
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // compute fv2
            real tFv2 = compute_fv2(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute dfv2du
            Matrix< DDRMat > tdfv2du;
            compute_dfv2du(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aDofTypes,
                    tdfv2du );

            // compute dsbardu
            adsbardu =
                    tFIViscosity->val() * tdfv2du /
                    std::pow( mKappa * tWallDistance, 2.0 );

            // if derivative dof type is viscosity dof type
            if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
            {
                // add contribution
                adsbardu += tFv2 * tFIViscosity->N() /
                        std::pow( mKappa * tWallDistance, 2.0 );
            }

            // if wall distance depends on derivative dof type
            if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                    ( tWallDistance > mWallDistanceEpsilon ) )
            {
                // add contribution to dsbardu
                adsbardu -=
                        2.0 * tFv2 * tFIViscosity->val()( 0 ) * aPropWallDistance->dPropdDOF( aDofTypes ) /
                        ( std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
            }
        }

        //------------------------------------------------------------------------------

        real compute_smod(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // compute s
            real tS = compute_s(
                    aVelocityDofGroup,
                    aMasterFIManager );

            // compute SBar
            real tSBar = compute_sbar(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute s
            real tSMod = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) /
                    ( ( mCv3 - 2.0 * mCv2 ) * tS - tSBar );

            return tSMod;
        }

        //------------------------------------------------------------------------------

        void compute_dsmoddu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsmoddu )
        {
            // compute s
            real tS = compute_s(
                    aVelocityDofGroup,
                    aMasterFIManager );

            // compute dsbardu
            Matrix< DDRMat > tdsdu;
            compute_dsdu(
                    aVelocityDofGroup,
                    aMasterFIManager ,
                    aDofTypes,
                    tdsdu );

            // compute SBar
            real tSBar = compute_sbar(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute dsdu
            Matrix< DDRMat > tdsbardu;
            compute_dsbardu(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance,
                    aDofTypes,
                    tdsbardu );

            // compute smod num
            real tSModNum = tS * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar );

            // compute smod deno
            real tSModDeno = ( mCv3 - 2.0 * mCv2 ) * tS - tSBar;

            // compute dsmoddu
            adsmoddu = ( ( tdsdu * ( std::pow( mCv2, 2 ) * tS + mCv3 * tSBar ) +
                    tS * ( std::pow( mCv2, 2 ) * tdsdu + mCv3 * tdsbardu ) ) * tSModDeno -
                    tSModNum * ( ( mCv3 - 2.0 * mCv2 ) * tdsdu - tdsbardu ) ) /
                            std::pow( tSModDeno, 2 );
        }

        //------------------------------------------------------------------------------

        real compute_stilde(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // compute s
            real tS = compute_s(
                    aVelocityDofGroup,
                    aMasterFIManager );

            // compute SBar
            real tSBar = compute_sbar(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // init stilde
            real tSTilde = tS;

            // compute STilde
            if( tSBar >= - mCv2 * tS )
            {
                tSTilde += tSBar;
            }
            else
            {
                // compute sMod
                tSTilde += compute_smod(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );
            }

            // return STilde
            return tSTilde;
        }

        //------------------------------------------------------------------------------

        void compute_dstildedu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adstildedu )
        {
            // compute S
            real tS = compute_s(
                    aVelocityDofGroup,
                    aMasterFIManager );

            // compute dSdu
            Matrix< DDRMat > tdSdu;
            compute_dsdu(
                    aVelocityDofGroup,
                    aMasterFIManager ,
                    aDofTypes,
                    tdSdu );

            // compute SBar
            real tSBar = compute_sbar(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // init dstildedu
            adstildedu = tdSdu;

            // compute dstildedu
            if( tSBar >= - mCv2 * tS )
            {
                // compute dSdu
                Matrix< DDRMat > tdSBardu;
                compute_dsbardu(
                        aViscosityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdSBardu );

                // add dSbardu
                adstildedu += tdSBardu;
            }
            else
            {
                // compute dSModdu
                Matrix< DDRMat > tdSModdu;
                compute_dsmoddu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdSModdu );

                // compute sMod
                adstildedu += tdSModdu;
            }
        }

        //------------------------------------------------------------------------------

        real compute_ft2(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity )
        {
            // compute chi
            real tChi = compute_chi(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

            // compute ft2
            real tFt2 = mCt3 * std::exp( - mCt4 * std::pow( tChi, 2.0 ) );

            return tFt2;
        }

        //------------------------------------------------------------------------------

        void compute_dft2du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adft2du )
        {
            // compute ft2
            real tFt2 = compute_ft2(
                    aViscosityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity );

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

            // compute dft2du
            adft2du = - mCt4 * tFt2 * 2.0 * tChi * tdchidu;
        }

        //------------------------------------------------------------------------------

        real compute_r(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // compute stilde
            real tSTilde = compute_stilde(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // threshold stilde
            tSTilde = std::max( tSTilde, mEpsilon );

            // get the residual dof FI (here viscosity)
            Field_Interpolator * tFIViscosity =
                    aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

            // get the wall distance value
            real tWallDistance = aPropWallDistance->val()( 0 );

            // threshold wall distance
            tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

            // compute viscosity / ( stilde * kappa² * d² )
            real tR = tFIViscosity->val()( 0 ) / ( tSTilde * std::pow( mKappa * tWallDistance, 2.0 ) );

//            // check that r is finite and greater than zero or set it to mRLim
//            Matrix<DDRMat> tRMatrix( 1, 1, tR );
//            Matrix<DDRMat> tInvRMatrix( 1, 1, 1/tR );
//            if( !isfinite( tRMatrix ) || !isfinite(tInvRMatrix) )
//            {
//                tR = mRLim;
//            }

            return std::min( tR, mRLim );
        }

        //------------------------------------------------------------------------------

        void compute_drdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adrdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adrdu
            adrdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients() );

            // compute r
            real tR = compute_r(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // if r < 10
            if( tR < mRLim )
            {
                // get the viscosity dof FI
                Field_Interpolator * tFIViscosity =
                        aMasterFIManager->get_field_interpolators_for_type( aViscosityDofGroup( 0 ) );

                // get the wall distance value
                real tWallDistance = aPropWallDistance->val()( 0 );

                // threshold wall distance
                tWallDistance = std::max( tWallDistance, mWallDistanceEpsilon );

                // compute stilde
                real tSTilde = compute_stilde(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance );

                // threshold stilde
                tSTilde = std::max( tSTilde, mEpsilon );

                // compute dStildedu
                Matrix< DDRMat > tdSTildedu;
                compute_dstildedu(
                        aViscosityDofGroup,
                        aVelocityDofGroup,
                        aMasterFIManager,
                        aPropKinViscosity,
                        aPropWallDistance,
                        aDofTypes,
                        tdSTildedu );

                // add contribution from dStildedu
                adrdu = - tFIViscosity->val() * tdSTildedu / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );

                // if dof type is viscosity
                if( aDofTypes( 0 ) == aViscosityDofGroup( 0 ) )
                {
                    // add contribution from viscosity
                    adrdu += tSTilde * tDerFI->N() / std::pow( tSTilde * mKappa * tWallDistance, 2.0 );
                }

                // if wall distance depends on derivative dof type
                if( ( aPropWallDistance->check_dof_dependency( aDofTypes ) ) &&
                        ( tWallDistance > mWallDistanceEpsilon ) )
                {
                    // add contribution from wall distance
                    adrdu -= 2.0 * tFIViscosity->val()( 0 ) * aPropWallDistance->dPropdDOF( aDofTypes ) /
                            ( tSTilde * std::pow( mKappa, 2.0 ) * std::pow( tWallDistance, 3.0 ) );
                }
            }
            else
            {
                adrdu.fill( 0.0 );
            }
        }

        //------------------------------------------------------------------------------

        real compute_g(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // compute r
            real tR = compute_r(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute g
            return tR + mCw2 * ( std::pow( tR, 6.0 ) - tR );
        }

        //------------------------------------------------------------------------------

        void compute_dgdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adgdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adgdu
            adgdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients() );

            // compute r
            real tR = compute_r(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute drdu
            Matrix< DDRMat > tdrdu;
            compute_drdu(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance,
                    aDofTypes,
                    tdrdu );

            // compute adgdu
            adgdu = ( 1.0 + mCw2 * ( 6.0 * std::pow( tR, 5.0 ) - 1.0 ) ) * tdrdu;
        }

        //------------------------------------------------------------------------------

        real compute_fw(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance )
        {
            // compute g
            real tG = compute_g(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute fw
            real tFw = ( 1.0 + std::pow( mCw3, 6.0 ) ) /
                    ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) );
            tFw = tG * std::pow( tFw, 1.0 / 6.0 );

            return tFw;
        }

        //------------------------------------------------------------------------------

        void compute_dfwdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfwdu )
        {
            // get the derivative dof FIs
            Field_Interpolator * tDerFI =
                    aMasterFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // init adfwdu
            adfwdu.set_size( 1, tDerFI->get_number_of_space_time_coefficients() );

            // compute g
            real tG = compute_g(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // compute dgdu
            Matrix< DDRMat > tdgdu;
            compute_dgdu(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance,
                    aDofTypes,
                    tdgdu );

            // compute fw
            real tFw = compute_fw(
                    aViscosityDofGroup,
                    aVelocityDofGroup,
                    aMasterFIManager,
                    aPropKinViscosity,
                    aPropWallDistance );

            // init adfwdu
            adfwdu = ( tFw * std::pow( mCw3, 6.0 ) * tdgdu ) /
                    ( tG * ( std::pow( tG, 6.0 ) + std::pow( mCw3, 6.0 ) ) );
        }

    } /* namespace fem */
} /* namespace moris */

