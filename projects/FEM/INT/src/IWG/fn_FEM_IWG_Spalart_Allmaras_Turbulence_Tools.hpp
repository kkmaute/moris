/*
 * fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp
 *
 *  Created on: Apr 15, 2021
 *      Author: noel
 * 
 * Free functions used in the IWGs for the Spalart-Allmaras turbulence model
 */

#ifndef SRC_FEM_FN_FEM_IWG_SPALART_ALLMARAS_TOOLS_HPP_
#define SRC_FEM_FN_FEM_IWG_SPALART_ALLMARAS_TOOLS_HPP_

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Property.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------
        /**
         * compute production term
         * P = cb1 * ( 1 - ft2 ) * STilde * modViscosity
         * if modViscosity >= 0
         * P = cb1 * ( 1 - ct3 ) * S * modViscosity
         * if modViscosity <  0
         */
        real compute_production_term(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >        & aPropKinViscosity,
                const std::shared_ptr< Property >        & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the production term
         * @param[ in ] aDofTypes      a list of dof type wrt which
         *                             the derivative is requested
         * @param[ in ] adproductiondu a matrix to fill with dproductiondu
         */
        void compute_dproductiontermdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adproductiondu );

        //------------------------------------------------------------------------------
        /**
         * compute production coefficient
         * P = cb1 * ( 1 - ft2 ) * STilde
         * if modViscosity >= 0
         * P = cb1 * ( 1 - ct3 ) * S
         * if modViscosity <  0
         */
        real compute_production_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the production coefficient
         * @param[ in ] aDofTypes      a list of dof type wrt which
         *                             the derivative is requested
         * @param[ in ] adproductiondu a matrix to fill with dproductiondu
         */
        void compute_dproductioncoefficientdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adproductiondu );

        //------------------------------------------------------------------------------
        /**
         * compute wall destruction term
         * D = ( cw1 * fw - cb1 * ft2 / kappa² ) * ( modViscosity / wallD )²
         * if modViscosity >= 0
         * D = - cw1 * ( modViscosity / wallD )²
         * if modViscosity <  0
         */
        real compute_wall_destruction_term(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the wall destruction term
         * @param[ in ] aDofTypes           a list of dof type wrt which
         *                                  the derivative is requested
         * @param[ in ] adwalldestructiondu a matrix to fill with dwalldestructiondu
         */
        void compute_dwalldestructiontermdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwalldestructiondu );


        //------------------------------------------------------------------------------
        /**
         * compute wall destruction coefficient
         * D = ( cw1 * fw - cb1 * ft2 / kappa² ) * modViscosity / wallD²
         * if modViscosity >= 0
         * D = - cw1 * modViscosity / wallD²
         * if modViscosity <  0
         */
        real compute_wall_destruction_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the wall destruction coefficient
         * @param[ in ] aDofTypes           a list of dof type wrt which
         *                                  the derivative is requested
         * @param[ in ] adwalldestructiondu a matrix to fill with dwalldestructiondu
         */
        void compute_dwalldestructioncoefficientdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adwalldestructiondu );

        //------------------------------------------------------------------------------
        /**
         * compute diffusion coefficient
         * Diff = kinViscosity + modViscosity
         * if modViscosity >= 0
         * Diff = kinViscosity + modViscosity * fn
         * if modViscosity <  0
         */
        real compute_diffusion_coefficient(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the diffusion coefficient
         * @param[ in ] aDofTypes     a list of dof type wrt which
         *                            the derivative is requested
         * @param[ in ] addiffusiondu a matrix to fill with ddiffusiondu
         */
        void compute_ddiffusiondu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondu );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the diffusion coefficient
         * @param[ in ] aDofTypes     a list of dof type wrt which
         *                            the derivative is requested
         * @param[ in ] addiffusiondx a matrix to fill with addiffusiondx
         */
        void compute_ddiffusiondx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & addiffusiondx );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the space derivative of diffusion coefficient
         * @param[ in ] aDofTypes     a list of dof type wrt which
         *                            the derivative is requested
         * @param[ in ] addiffusiondxdu a matrix to fill with ddiffusiondxdu
         */
        void compute_ddiffusiondxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & addiffusiondxdu );

        //------------------------------------------------------------------------------
        /**
         * compute fn
         * fn = ( cn1 + chi³ ) / ( cn1 - chi³)
         */
        real compute_fn(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of fn wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adfndu     a matrix to fill with dfndu
         */
        void compute_dfndu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfndu );

        //------------------------------------------------------------------------------
        /**
         * compute the space derivative of fn
         * @param[ in ] adfndx     a matrix to fill with dfndx
         */
        void compute_dfndx(
                        const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                        Field_Interpolator_Manager         * aMasterFIManager,
                        const std::shared_ptr< Property >  & aPropKinViscosity,
                        Matrix< DDRMat >                   & adfndx );
        //------------------------------------------------------------------------------
             /**
              * compute the derivative of space derivative of fn wrt to a dof type
              * @param[ in ] aDofTypes  a list of dof type wrt which
              *                         the derivative is requested
              * @param[ in ] adfndxdu   a matrix to fill with dfndxdu
              */
        void compute_dfndxdu(
                     const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                     Field_Interpolator_Manager         * aMasterFIManager,
                     const std::shared_ptr< Property >  & aPropKinViscosity,
                     const moris::Cell< MSI::Dof_Type > & aDofTypes,
                     Matrix< DDRMat >                   & adfndxdu );

        //------------------------------------------------------------------------------
        /**
         * compute Wij = 0.5 * ( dui/dxj - duj/dxi )
         * @param[ out ] Wij
         */
        void compute_wij(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                Matrix< DDRMat >                   & aWij );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of Wij wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adwijdu    a matrix to fill with dwijdu
         */
        void compute_dwijdu(
                const moris::Cell< MSI::Dof_Type >  & aVelocityDofGroup,
                Field_Interpolator_Manager          * aMasterFIManager,
                const  moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                    & adwijdu );

        //------------------------------------------------------------------------------
        /**
         * compute chi = viscosityDof / viscosityPtop
         * @param[ out ] chi
         */
        real compute_chi(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of chi wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adchidu    a matrix to fill with dchidu
         */
        void compute_dchidu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidu );

        //------------------------------------------------------------------------------
        /**
         * compute the space derivative of chi
         * @param[ in ] adchidu    a matrix to fill with dchidu
         */
        void compute_dchidx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & adchidx );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of the space derivative of chi wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adchidxdu    a matrix to fill with dchidu
         */
        void compute_dchidxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adchidxdu );

        //------------------------------------------------------------------------------
        /**
         * compute fv1 = chi³ / ( chi³ + cv1³)
         * @param[ out ] fv1
         */
        real compute_fv1(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of fv1 wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adfv1du    a matrix to fill with dfv1du
         */
        void compute_dfv1du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1du );

        //------------------------------------------------------------------------------
        /**
         * compute space derivative of fv1 = chi³ / ( chi³ + cv1³)
         * @param[ out ] adfv1dx a matrix to fill with adfv1dx
         */
        void compute_dfv1dx(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                Matrix< DDRMat >                   & adfv1dx );

        //------------------------------------------------------------------------------
             /**
              * compute the derivative of the space derivative of fv1 wrt to a dof type
              * @param[ in ] aDofTypes  a list of dof type wrt which
              *                         the derivative is requested
              * @param[ in ] adfv1dxdu  a matrix to fill with adfv1dxdu
              */
        void compute_dfv1dxdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv1dxdu );

        //------------------------------------------------------------------------------
        /**
         * compute fv2 = 1 - chi / ( 1 + chi * fv1 )
         * @param[ out ] fv2
         */
        real compute_fv2(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of fv2 wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adfv2du    a matrix to fill with dfv2du
         */
        void compute_dfv2du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfv2du );

        //------------------------------------------------------------------------------
        /**
         * compute SBar = sqrt( 2 * Wij * Wij )
         * @param[ out ] SBar
         */
        real compute_sbar(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );
        //------------------------------------------------------------------------------
        /**
         * compute the derivative of SBar wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adSBardu a matrix to fill with dSBardu
         */
        void compute_dsbardu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsbardu );

        //------------------------------------------------------------------------------
        /**
         * compute S = viscosity * fv2 / ( kappa² * d² )
         * @param[ out ] S
         */
        real compute_s(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of S wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adSdu a matrix to fill with dSdu
         */
        void compute_dsdu(
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsdu );

        //------------------------------------------------------------------------------
        /**
         * compute Smod =
         * @param[ out ] Smod
         */
        real compute_smod(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of Smod wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adSmoddu a matrix to fill with dSmoddu
         */
        void compute_dsmoddu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adsmoddu );

        //------------------------------------------------------------------------------
        /**
         * compute STilde = SBar + viscosity * fv2 / ( kappa² * d² )
         * @param[ out ] STilde
         */
        real compute_stilde(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of STilde wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adSTildedu a matrix to fill with dSTildedu
         */
        void compute_dstildedu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adstildedu );

        //------------------------------------------------------------------------------
        /**
         * compute ft2 = ct3 * exp( -ct4 * chi²), chi = viscosityDof/viscosityProp
         * @param[ out ] ft2
         */
        real compute_ft2(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of ft2 wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adft2du a matrix to fill with dft2du
         */
        void compute_dft2du(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adft2du );

        //------------------------------------------------------------------------------
        /**
         * compute r = min( viscosity / ( STilde * kappa² * d² ), 10 )
         * @param[ out ] r
         */
        real compute_r(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of r wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adrdu a matrix to fill with drdu
         */
        void compute_drdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adrdu );

        //------------------------------------------------------------------------------
        /**
         * compute g = r + cw2 * ( r⁶ - r )
         * @param[ out ] g
         */
        real compute_g(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of g wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adgdu a matrix to fill with dgdu
         */
        void compute_dgdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adgdu );

        //------------------------------------------------------------------------------
        /**
         * compute fw = g * ( ( 1 + cw3⁶ ) / ( g⁶ + cw3⁶ ) )^(1/6)
         * g = r + cw2 * ( r⁶ - r )
         * r = min( viscosity / ( STilde * kappa² * d² ), 10 )
         * @param[ out ] fw
         */
        real compute_fw(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance );

        //------------------------------------------------------------------------------
        /**
         * compute the derivative of fw wrt to a dof type
         * @param[ in ] aDofTypes  a list of dof type wrt which
         *                         the derivative is requested
         * @param[ in ] adfwdu a matrix to fill with dfwdu
         */
        void compute_dfwdu(
                const moris::Cell< MSI::Dof_Type > & aViscosityDofGroup,
                const moris::Cell< MSI::Dof_Type > & aVelocityDofGroup,
                Field_Interpolator_Manager         * aMasterFIManager,
                const std::shared_ptr< Property >  & aPropKinViscosity,
                const std::shared_ptr< Property >  & aPropWallDistance,
                const moris::Cell< MSI::Dof_Type > & aDofTypes,
                Matrix< DDRMat >                   & adfwdu );


    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_IWG_SPALART_ALLMARAS_TOOLS_HPP_ */
