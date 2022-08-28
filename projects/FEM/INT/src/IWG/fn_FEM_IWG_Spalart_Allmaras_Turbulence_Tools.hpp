/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.hpp
 *
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

    } /* namespace fem */
} /* namespace moris */

#endif /* SRC_FEM_FN_FEM_IWG_SPALART_ALLMARAS_TOOLS_HPP_ */

