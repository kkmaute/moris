/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Spalart_Allmaras_Turbulence_Tools.cpp
 *
 */

#include "fn_FEM_IWG_Crosswind_Stabilization_Tools.hpp"
#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "fn_eye.hpp"
#include "fn_clip_value.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void compute_cgradxw(
            const moris::Vector< MSI::Dof_Type > & aVelocityDofGroup,
            const moris::Vector< MSI::Dof_Type > & aTargetDofGroup,
            Field_Interpolator_Manager         * aLeaderFIManager,
            const real                         & aSpaceDim,
            const real                         & aEpsilon,
            bool                               & aIsCrosswind,
            Matrix< DDRMat >                   & acgradxw )
        {
            // get the target dof FI
            Field_Interpolator* tFITarget =
                    aLeaderFIManager->get_field_interpolators_for_type( aTargetDofGroup( 0 ) );

            // if crosswind
            if( aIsCrosswind )
            {
                // get the velocity dof FI
                Field_Interpolator* tFIVelocity = //
                        aLeaderFIManager->get_field_interpolators_for_type( aVelocityDofGroup( 0 ) );

                // compute the square of the norm of the velocity
                const real tNormVelocitySquare = //
                        std::max( dot( tFIVelocity->val(), tFIVelocity->val() ), aEpsilon );

                // build an identity matrix
                eye( aSpaceDim, aSpaceDim, acgradxw );

                // add crosswind
                acgradxw -= tFIVelocity->val() * tFIVelocity->val_trans() / tNormVelocitySquare;

                // compute c * gradx(w)
                acgradxw = acgradxw * tFITarget->gradx( 1 );

                // vectorize c * gradx(w)
                acgradxw = vectorize( acgradxw );

            }
            // if isotropic diffusion
            else
            {
                // compute I * gradx(w)
                acgradxw = vectorize( tFITarget->gradx( 1 ) );
            }
        }

        //------------------------------------------------------------------------------

        void compute_dcgradxwdu(
                const moris::Vector< MSI::Dof_Type > & aVelocityDofGroup,
                const moris::Vector< MSI::Dof_Type > & aTargetDofGroup,
                Field_Interpolator_Manager         * aLeaderFIManager,
                const moris::Vector< MSI::Dof_Type > & aDofTypes,
                const real                         & aSpaceDim,
                const real                         & aEpsilon,
                bool                               & aIsCrosswind,
                Matrix< DDRMat >                   & adcgradxwdu )
        {
            // get the dof type FI
            Field_Interpolator* tFIDer =
                    aLeaderFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

            // get the number of coefficients
            uint tNumCoeff = tFIDer->get_number_of_space_time_coefficients();

            // get the FI for scalar field
            Field_Interpolator* tFITarget =
                    aLeaderFIManager->get_field_interpolators_for_type( aTargetDofGroup( 0 ) );

            // get number of coefficient per field
            uint tTargetDim = tFITarget->get_number_of_fields();
            uint tNumCoeffPerField = tNumCoeff / tTargetDim;

            // init adcgradxwdu
            adcgradxwdu.set_size( aSpaceDim * tTargetDim, tNumCoeff, 0.0 );

            // if crosswind diffuusion
            if( aIsCrosswind )
            {
                // get the FI for scalar field
                Field_Interpolator* tFIVelocity = //
                        aLeaderFIManager->get_field_interpolators_for_type( aVelocityDofGroup( 0 ) );

                // compute the square of the norm of the velocity
                const real tNormVelocitySquare = //
                        std::max( dot( tFIVelocity->val(), tFIVelocity->val() ), aEpsilon );

                // if derivative dof is velocity
                if( aDofTypes( 0 ) == aVelocityDofGroup( 0 ) )
                {
                    // compute the derivative of the square of the norm of the velocity
                    Matrix< DDRMat > tdNormVelocitySquaredu( 1, tNumCoeff );
                    if( tNormVelocitySquare > aEpsilon )
                    {
                        tdNormVelocitySquaredu = - 2.0 * tFIVelocity->val_trans() * tFIVelocity->N() / std::pow( tNormVelocitySquare, 2.0 );
                    }
                    else
                    {
                        tdNormVelocitySquaredu.fill( 0.0 );
                    }

                    // compute the derivative of the crosswind orthogonal projection
                    Matrix< DDRMat > tdCdu( aSpaceDim * aSpaceDim, tNumCoeff );

                    // build derivatives of component of the crosswind projection
                    for( uint i = 0; i < aSpaceDim; i++ )
                    {
                        for( uint j = 0; j < aSpaceDim; j++ )
                        {
                            // get the row index
                            uint tRowIndex = i * aSpaceDim + j;

                            // fill contribution
                            tdCdu( { tRowIndex, tRowIndex }, { 0, tNumCoeff - 1 } ) = //
                                    - tFIVelocity->N()({ i, i }, { 0, tNumCoeff - 1 } ) * tFIVelocity->val()( j ) / tNormVelocitySquare //
                                    - tFIVelocity->val()( i ) * tFIVelocity->N()({ j, j }, { 0, tNumCoeff - 1 } ) / tNormVelocitySquare //
                                    - tFIVelocity->val()( i ) * tFIVelocity->val()( j ) * tdNormVelocitySquaredu;
                        }

                        // loop over the number of fields
                        for( uint k = 0; k < tTargetDim; k++ )
                        {
                            // get the row index
                            uint tRowIndex = k * tTargetDim + i;

                            // add contribution to derivative
                            adcgradxwdu( //
                                    { tRowIndex, tRowIndex }, //
                                    { 0, tNumCoeff - 1 } ) += //
                                            trans( tFITarget->gradx( 1 )( { 0, aSpaceDim - 1 }, { k, k } ) ) * //
                                            tdCdu( { i * aSpaceDim, ( i + 1 ) * aSpaceDim - 1 }, { 0, tNumCoeff - 1 } );
                        }
                    }
                }

                // if derivative dof is target
                if( aDofTypes( 0 ) == aTargetDofGroup( 0 ) )
                {
                    // build crosswind orthogonal projection
                    Matrix< DDRMat > tC;

                    // build an identity matrix
                    eye( aSpaceDim, aSpaceDim, tC );

                    // build crosswind term
                    tC -= tFIVelocity->val() * tFIVelocity->val_trans() / tNormVelocitySquare;

                    // build derivatives of component of the crosswind projection
                    for( uint i = 0; i < aSpaceDim; i++ )
                    {
                        // loop over the number of field
                        for( uint k = 0; k < tTargetDim; k++ )
                        {
                            // get the row index
                            uint tRowIndex = k * aSpaceDim + i;

                            // add contribution to derivative
                            adcgradxwdu( //
                                    { tRowIndex, tRowIndex }, //
                                    { k * tNumCoeffPerField, ( k + 1 ) * tNumCoeffPerField - 1 } ) += //
                                            tC( { i, i }, { 0, aSpaceDim - 1 } ) * tFIDer->dnNdxn( 1 );
                        }
                    }
                }
            }
            // if isotropic diffusion
            else
            {
                // if derivative dof is target
                if( aDofTypes( 0 ) == aTargetDofGroup( 0 ) )
                {
                    // loop over the number of field
                    for( uint k = 0; k < tTargetDim; k++ )
                    {
                        // get the row index
                        uint tRowStartIndex = k * aSpaceDim;
                        uint tRowStopIndex  = ( k + 1 ) * aSpaceDim - 1;

                        // add contribution to derivative
                        adcgradxwdu(                                                              //
                                { tRowStartIndex, tRowStopIndex },                                         //
                                { k * tNumCoeffPerField, ( k + 1 ) * tNumCoeffPerField - 1 } ) += //
                                        tFIDer->dnNdxn( 1 ).matrix_data();
                    }
                }
            }
        }

    } /* namespace fem */
} /* namespace moris */

