/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_IWG_Compressible_NS_Flux_Matrices_Dof_Derivs.cpp
 *
 */

#include "fn_FEM_IWG_Compressible_NS.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    void eval_A0_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA0dDOF )
    {
        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::eval_A0_DOF - list of aResidualDofTypes not supported, see messages above." );
        MORIS_ASSERT( ( aMM != nullptr ) and ( aCM != nullptr ) and ( aLeaderFIManager != nullptr ),
                "fn_FEM_IWG_Compressible_NS::eval_A0_DOF - nullptr provided in inputs check MM, CM, and FI Manager provided." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // get commonly used values
        real tRho    = aMM->density()( 0 );
        real tAlphaP = aMM->AlphaP()( 0 );
        real tBetaT  = aMM->BetaT()( 0 );
        real tCv     = aMM->Cv()( 0 );
        real tEtot   = aCM->Energy()( 0 );
        real tUx     = tFIVelocity->val()( 0 );
        real tUy     = tFIVelocity->val()( 1 );
        real tUz     = 0.0;
        if ( tNumSpaceDims == 3 )
        {
            tUz = tFIVelocity->val()( 2 );
        }

        // divide the N-vector for the velocity
        Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
        Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
        Matrix< DDRMat > tNUz;
        if ( tNumSpaceDims == 3 )
        {
            tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
        }

        // variable index for the third variable depending on spatial dimension
        uint tThirdVarIndex = tNumSpaceDims + 1;

        // =======================
        // Assemble A0 derivatives
        // =======================
        // clang-format off

            // derivative matrix for FIRST ROW OF A0

            adA0dDOF( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) );

            adA0dDOF( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) );

            adA0dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) );

            adA0dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) );

            // derivative matrix for SECOND ROW OF A0

            adA0dDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA0dDOF( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tBetaT * tNUx;

            adA0dDOF( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA0dDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

            adA0dDOF( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

            adA0dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tUx * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA0dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tNUx;

            adA0dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for THIRD ROW OF A0

            adA0dDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA0dDOF( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tBetaT * tNUy;

            adA0dDOF( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA0dDOF( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

            adA0dDOF( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

            adA0dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tUy * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA0dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tNUy;

            adA0dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for FOURTH ROW OF A0 IF 3D
            if ( tNumSpaceDims == 3 )
            {
                adA0dDOF( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

                adA0dDOF( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        tRho * tBetaT * tNUz;

                adA0dDOF( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

                adA0dDOF( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                        aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

                adA0dDOF( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

                adA0dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                        tUz * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

                adA0dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        -1.0 * tRho * tAlphaP * tNUz;

                adA0dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUz * ( -1.0 * tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) - tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );
            }

            // derivative matrix for LAST ROW OF A0

            adA0dDOF( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tEtot * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 0 ) );

            adA0dDOF( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) =
                    tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA0dDOF( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tEtot * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 2 ) );

            adA0dDOF( tThirdVarIndex )( { 1, tNumSpaceDims }, { 0, tNumBases - 1 } ) =
                    tFIVelocity->val() * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA0dDOF( tThirdVarIndex )( { 1, tNumSpaceDims }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) =
                    tRho * tFIVelocity->N();

            adA0dDOF( tThirdVarIndex )( { 1, tNumSpaceDims }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tFIVelocity->val() * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA0dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tCv * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 0 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 0 ) );

            adA0dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) =
                    -1.0 * tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA0dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tCv * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 2 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 2 ) );

        // clang-format on
    }

    //------------------------------------------------------------------------------

    void eval_A1_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA1dDOF )
    {
        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::eval_A1_DOF - list of aResidualDofTypes not supported, see messages above." );
        MORIS_ASSERT( ( aMM != nullptr ) and ( aCM != nullptr ) and ( aLeaderFIManager != nullptr ),
                "fn_FEM_IWG_Compressible_NS::eval_A1_DOF - nullptr provided in inputs check MM, CM, and FI Manager provided." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // get commonly used values
        real tRho    = aMM->density()( 0 );
        real tAlphaP = aMM->AlphaP()( 0 );
        real tBetaT  = aMM->BetaT()( 0 );
        real tCv     = aMM->Cv()( 0 );
        real tEtot   = aCM->Energy()( 0 );
        real tUx     = tFIVelocity->val()( 0 );
        real tUy     = tFIVelocity->val()( 1 );
        real tUz     = 0.0;
        if ( tNumSpaceDims == 3 )
        {
            tUz = tFIVelocity->val()( 2 );
        }

        // divide the N-vector for the velocity
        Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
        Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
        Matrix< DDRMat > tNUz;
        if ( tNumSpaceDims == 3 )
        {
            tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
        }

        // variable index for the third variable depending on spatial dimension
        uint tThirdVarIndex = tNumSpaceDims + 1;

        // =======================
        // Assemble A1 derivatives
        // =======================
        // clang-format off

            // derivative matrix for FIRST ROW OF A1

            adA1dDOF( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 0 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tBetaT * tNUx;

            adA1dDOF( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA1dDOF( 0 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

            adA1dDOF( 0 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

            adA1dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUx * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tNUx;

            adA1dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUx * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for SECOND ROW OF A1

            adA1dDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    2.0 * tRho * tBetaT * tUx * tNUx;

            adA1dDOF( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * tUx * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA1dDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    2.0 * tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA1dDOF( 1 )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    2.0 * tRho * tNUx;

            adA1dDOF( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    2.0 * tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA1dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUx * tUx * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -2.0 * tRho * tAlphaP * tUx * tNUx;

            adA1dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUx * tUx * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for THIRD ROW OF A1

            adA1dDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 2 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUy * tNUx;

            adA1dDOF( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUx * tNUy;

            adA1dDOF( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA1dDOF( 2 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA1dDOF( 2 )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tNUy;

            adA1dDOF( 2 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA1dDOF( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA1dDOF( 2 )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tNUx;

            adA1dDOF( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA1dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUx * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUy * tNUx;

            adA1dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUx * tNUy;

            adA1dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUx * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for FOURTH ROW OF A1 IF 3D
            if ( tNumSpaceDims == 3 )
            {
                adA1dDOF( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tUx * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

                adA1dDOF( 3 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                        tRho * tBetaT * tUz * tNUx;

                adA1dDOF( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        tRho * tBetaT * tUx * tNUz;

                adA1dDOF( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUx * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

                adA1dDOF( 3 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA1dDOF( 3 )( { 1, 1 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        tRho * tNUz;

                adA1dDOF( 3 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

                adA1dDOF( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                        tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA1dDOF( 3 )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) =
                        tRho * tNUx;

                adA1dDOF( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

                adA1dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                        -1.0 * tUx * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

                adA1dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                        -1.0 * tRho * tAlphaP * tUz * tNUx;

                adA1dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        -1.0 * tRho * tAlphaP * tUx * tNUz;

                adA1dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        -1.0 * tUx * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );
            }

            // derivative matrix for LAST ROW OF A1

            adA1dDOF( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) );

            adA1dDOF( tThirdVarIndex )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    ( tBetaT * tEtot + 1.0 ) * tNUx;

            adA1dDOF( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    tUx * tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA1dDOF( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

            adA1dDOF( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) + aMM->PressureDOF( aResidualDofTypes( 0 ) ) + tUx * tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA1dDOF( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =
                   2.0 * tRho * tUx * tNUx;

            adA1dDOF( tThirdVarIndex )( { 1, 1 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    aCM->dEnergydDOF( aResidualDofTypes( 1 ) ).matrix_data();

            adA1dDOF( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) + aMM->PressureDOF( aResidualDofTypes( 2 ) ) + tUx * tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA1dDOF( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                   tUx * tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA1dDOF( tThirdVarIndex )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =
                   tRho * tUy * tNUx;

            adA1dDOF( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                   tRho * tUx * tNUy;

            adA1dDOF( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   tUx * tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            if ( tNumSpaceDims == 3 ) // third velocity component only for 3D
            {
                adA1dDOF( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    tUx * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA1dDOF( tThirdVarIndex )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tUz * tNUx;

                adA1dDOF( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tUx * tNUz;

                adA1dDOF( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );
            }

            adA1dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tUx * ( tCv * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 0 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) ) ;

            adA1dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUx;

            adA1dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -=
                    tUx * tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA1dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * ( tCv * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 2 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

        // clang-format on
    }

    //------------------------------------------------------------------------------

    void eval_A2_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA2dDOF )
    {
        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::eval_A2_DOF - list of aResidualDofTypes not supported, see messages above." );
        MORIS_ASSERT( ( aMM != nullptr ) and ( aCM != nullptr ) and ( aLeaderFIManager != nullptr ),
                "fn_FEM_IWG_Compressible_NS::eval_A2_DOF - nullptr provided in inputs check MM, CM, and FI Manager provided." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // get commonly used values
        real tRho    = aMM->density()( 0 );
        real tAlphaP = aMM->AlphaP()( 0 );
        real tBetaT  = aMM->BetaT()( 0 );
        real tCv     = aMM->Cv()( 0 );
        real tEtot   = aCM->Energy()( 0 );
        real tUx     = tFIVelocity->val()( 0 );
        real tUy     = tFIVelocity->val()( 1 );
        real tUz     = 0.0;
        if ( tNumSpaceDims == 3 )
        {
            tUz = tFIVelocity->val()( 2 );
        }

        // divide the N-vector for the velocity
        Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
        Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
        Matrix< DDRMat > tNUz;
        if ( tNumSpaceDims == 3 )
        {
            tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );
        }

        // variable index for the third variable depending on spatial dimension
        uint tThirdVarIndex = tNumSpaceDims + 1;

        // =======================
        // Assemble A2 derivatives
        // =======================
        // clang-format off

            // derivative matrix for FIRST ROW OF A2

            adA2dDOF( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 0 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tBetaT * tNUy;

            adA2dDOF( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA2dDOF( 0 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

            adA2dDOF( 0 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

            adA2dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tNUy;

            adA2dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for SECOND ROW OF A2 (is equivalent ot third row of A1)

            adA2dDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUy * tNUx;

            adA2dDOF( 1 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUx * tNUy;

            adA2dDOF( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA2dDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA2dDOF( 1 )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tNUy;

            adA2dDOF( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA2dDOF( 1 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA2dDOF( 1 )( { 2, 2 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tNUx;

            adA2dDOF( 1 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA2dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUx * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUy * tNUx;

            adA2dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUx * tNUy;

            adA2dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUx * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for THIRD ROW OF A2

            adA2dDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUy * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    2.0 * tRho * tBetaT * tUy * tNUy;

            adA2dDOF( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * tUy * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA2dDOF( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    2.0 * tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA2dDOF( 2 )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    2.0 * tRho * tNUy;

            adA2dDOF( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    2.0 * tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA2dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUy * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -2.0 * tRho * tAlphaP * tUy * tNUy;

            adA2dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUy * tUy * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for FOURTH ROW OF A2
            if ( tNumSpaceDims == 3 )
            {
                adA2dDOF( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tUy * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

                adA2dDOF( 3 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                        tRho * tBetaT * tUz * tNUy;

                adA2dDOF( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        tRho * tBetaT * tUy * tNUz;

                adA2dDOF( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUy * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

                adA2dDOF( 3 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                        tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA2dDOF( 3 )( { 2, 2 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        tRho * tNUz;

                adA2dDOF( 3 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

                adA2dDOF( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                        tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA2dDOF( 3 )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                        tRho * tNUy;

                adA2dDOF( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

                adA2dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                        -1.0 * tUy * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

                adA2dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                        -1.0 * tRho * tAlphaP * tUz * tNUy;

                adA2dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                        -1.0 * tRho * tAlphaP * tUy * tNUz;

                adA2dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                        -1.0 * tUy * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );
            }

            // derivative matrix for LAST ROW OF A2

            adA2dDOF( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUy * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) );

            adA2dDOF( tThirdVarIndex )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    ( tBetaT * tEtot + 1.0 ) * tNUy;

            adA2dDOF( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    tUy * tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA2dDOF( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

            adA2dDOF( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                   tUx * tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA2dDOF( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =
                   tRho * tUy * tNUx;

            adA2dDOF( tThirdVarIndex )( { 1, 1 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                   tRho * tUx * tNUy;

            adA2dDOF( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   tUx * tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA2dDOF( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) + aMM->PressureDOF( aResidualDofTypes( 0 ) ) + tUy * tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA2dDOF( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                   2.0 * tRho * tUy * tNUy;

            adA2dDOF( tThirdVarIndex )( { 2, 2 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    aCM->dEnergydDOF( aResidualDofTypes( 1 ) ).matrix_data();

            adA2dDOF( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) + aMM->PressureDOF( aResidualDofTypes( 2 ) ) + tUy * tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            if ( tNumSpaceDims == 3 ) // third velocity component only for 3D
            {
                adA2dDOF( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    tUy * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

                adA2dDOF( tThirdVarIndex )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tUz * tNUy;

                adA2dDOF( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tUy * tNUz;

                adA2dDOF( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );
            }

            adA2dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tUy * ( tCv * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 0 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) ) ;

            adA2dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUy;

            adA2dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -=
                    tUy * tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA2dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * ( tCv * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 2 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

        // clang-format on
    }

    //------------------------------------------------------------------------------

    void eval_A3_DOF(
            const std::shared_ptr< Material_Model >     &aMM,
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            Vector< Matrix< DDRMat > >                  &adA3dDOF )
    {
        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::eval_A3_DOF - list of aResidualDofTypes not supported, see messages above." );
        MORIS_ASSERT( ( aMM != nullptr ) and ( aCM != nullptr ) and ( aLeaderFIManager != nullptr ),
                "fn_FEM_IWG_Compressible_NS::eval_A3_DOF - nullptr provided in inputs check MM, CM, and FI Manager provided." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // check number of space dimensions
        MORIS_ASSERT( tNumSpaceDims == 3, "fn_FEM_IWG_Compressible_NS::eval_A3_DOF - A3 only needed and defined for 3D." );

        // get commonly used values
        real tRho    = aMM->density()( 0 );
        real tAlphaP = aMM->AlphaP()( 0 );
        real tBetaT  = aMM->BetaT()( 0 );
        real tCv     = aMM->Cv()( 0 );
        real tEtot   = aCM->Energy()( 0 );
        real tUx     = tFIVelocity->val()( 0 );
        real tUy     = tFIVelocity->val()( 1 );
        real tUz     = tFIVelocity->val()( 2 );

        // divide the N-vector for the velocity
        Matrix< DDRMat > tNUx = tFIVelocity->N()( { 0, 0 }, { 0, tNumBases - 1 } );
        Matrix< DDRMat > tNUy = tFIVelocity->N()( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } );
        Matrix< DDRMat > tNUz = tFIVelocity->N()( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } );

        // variable index for the third variable depending on spatial dimension
        uint tThirdVarIndex = tNumSpaceDims + 1;

        // =======================
        // Assemble A3 derivatives
        // =======================
        // clang-format off

            // derivative matrix for FIRST ROW OF A3

            adA3dDOF( 0 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 0 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tBetaT * tNUz;

            adA3dDOF( 0 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA3dDOF( 0 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 0 ) ).matrix_data();

            adA3dDOF( 0 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    aMM->DensityDOF( aResidualDofTypes( 2 ) ).matrix_data();

            adA3dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tNUz;

            adA3dDOF( 0 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for SECOND ROW OF A3 (identical to fourth row of A1)

            adA3dDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUx * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 1 )( { 0, 0 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUz * tNUx;

            adA3dDOF( 1 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUx * tNUz;

            adA3dDOF( 1 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA3dDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                    tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( 1 )( { 1, 1 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tNUz;

            adA3dDOF( 1 )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( 1 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( 1 )( { 3, 3 }, { tNumBases, 2 * tNumBases - 1 } ) =
                    tRho * tNUx;

            adA3dDOF( 1 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUx * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUx * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, 2 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUz * tNUx;

            adA3dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUx * tNUz;

            adA3dDOF( 1 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUx * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for THIRD ROW OF A3 (identical to fourth row of A2)

            adA3dDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUy * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 2 )( { 0, 0 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUz * tNUy;

            adA3dDOF( 2 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tBetaT * tUy * tNUz;

            adA3dDOF( 2 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA3dDOF( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                    tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( 2 )( { 2, 2 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    tRho * tNUz;

            adA3dDOF( 2 )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( 2 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( 2 )( { 3, 3 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    tRho * tNUy;

            adA3dDOF( 2 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUy * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUy * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUz * tNUy;

            adA3dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    -1.0 * tRho * tAlphaP * tUy * tNUz;

            adA3dDOF( 2 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUy * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for FOURTH ROW OF A3

            adA3dDOF( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUz * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 3 )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    2.0 * tRho * tBetaT * tUz * tNUz;

            adA3dDOF( 3 )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * tUz * ( tBetaT * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) );

            adA3dDOF( 3 )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                    2.0 * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( 3 )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    2.0 * tRho * tNUz;

            adA3dDOF( 3 )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    2.0 * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    -1.0 * tUz * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    -2.0 * tRho * tAlphaP * tUz * tNUz;

            adA3dDOF( 3 )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    -1.0 * tUz * tUz * ( tAlphaP * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) );

            // derivative matrix for LAST ROW OF A3

            adA3dDOF( tThirdVarIndex )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                    tUz * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 0 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) );

            adA3dDOF( tThirdVarIndex )( { 0, 0 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    ( tBetaT * tEtot + 1.0 ) * tNUz;

            adA3dDOF( tThirdVarIndex )( { 0, 0 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    tUz * tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA3dDOF( tThirdVarIndex )( { 0, 0 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * ( tEtot * aMM->BetaTDOF( aResidualDofTypes( 2 ) ) + tBetaT * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

            adA3dDOF( tThirdVarIndex )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                   tUx * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( tThirdVarIndex )( { 1, 1 }, { tNumBases, 2 * tNumBases - 1 } ) =
                   tRho * tUz * tNUx;

            adA3dDOF( tThirdVarIndex )( { 1, 1 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                   tRho * tUx * tNUz;

            adA3dDOF( tThirdVarIndex )( { 1, 1 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   tUx * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( tThirdVarIndex )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                tUy * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( tThirdVarIndex )( { 2, 2 }, { 2 * tNumBases, 3 * tNumBases - 1 } ) =
                tRho * tUz * tNUy;

            adA3dDOF( tThirdVarIndex )( { 2, 2 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                tRho * tUy * tNUz;

            adA3dDOF( tThirdVarIndex )( { 2, 2 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                tUy * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( tThirdVarIndex )( { 3, 3 }, { 0, tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) + aMM->PressureDOF( aResidualDofTypes( 0 ) ) + tUz * tUz * aMM->DensityDOF( aResidualDofTypes( 0 ) );

            adA3dDOF( tThirdVarIndex )( { 3, 3 }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                   2.0 * tRho * tUz * tNUz;

            adA3dDOF( tThirdVarIndex )( { 3, 3 }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) +=
                    aCM->dEnergydDOF( aResidualDofTypes( 1 ) ).matrix_data();

            adA3dDOF( tThirdVarIndex )( { 3, 3 }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                   aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) + aMM->PressureDOF( aResidualDofTypes( 2 ) ) + tUz * tUz * aMM->DensityDOF( aResidualDofTypes( 2 ) );

            adA3dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 0, tNumBases - 1 } ) =
                    tUz * ( tCv * aMM->DensityDOF( aResidualDofTypes( 0 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 0 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 0 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 0 ) ) ) ;

            adA3dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { 3 * tNumBases, 4 * tNumBases - 1 } ) =
                    ( tRho * tCv - tAlphaP * tEtot ) * tNUz;

            adA3dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { tNumBases, ( tNumSpaceDims + 1 ) * tNumBases - 1 } ) -=
                    tUz * tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 1 ) );

            adA3dDOF( tThirdVarIndex )( { tThirdVarIndex, tThirdVarIndex }, { ( tNumSpaceDims + 1 ) * tNumBases, ( tNumSpaceDims + 2 ) * tNumBases - 1 } ) =
                    tUz * ( tCv * aMM->DensityDOF( aResidualDofTypes( 2 ) ) + tRho * aMM->CpDOF( aResidualDofTypes( 2 ) ) -
                    tEtot * aMM->AlphaPDOF( aResidualDofTypes( 2 ) ) - tAlphaP * aCM->dEnergydDOF( aResidualDofTypes( 2 ) ) );

        // clang-format on
    }

    //------------------------------------------------------------------------------

    void eval_KijYjDOF(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Vector< MSI::Dof_Type >               &aDofType,
            Vector< Matrix< DDRMat > >                  &aKijYjDOF )
    {
        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::eval_KijYjDOF_matrices - list of aResidualDofTypes not supported, see error messages above." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();
        if ( aDofType( 0 ) == MSI::Dof_Type::VX )
        {
            tNumBases = tNumBases * tNumSpaceDims;
        }

        // initialize dKijYj/dDOF, each cell entry represents the dof derivs of one column of Kij*Y,j
        aKijYjDOF.resize( tNumSpaceDims + 2 );

        // density / pressure residual - simply zero matrix
        aKijYjDOF( 0 ).set_size( tNumSpaceDims, tNumBases, 0.0 );

        // velocity residual
        // get the viscous stress Dof Deriv
        Matrix< DDRMat > tdTaudDOF = aCM->dFluxdDOF( aDofType, CM_Function_Type::MECHANICAL );

        // check that sizes match
        MORIS_ASSERT( ( tdTaudDOF.n_rows() == 3 * tNumSpaceDims - 3 ) and ( tdTaudDOF.n_cols() == tNumBases ),
                "fn_FEM_IWG_Compressible_NS::eval_KijYjDOF_matrices - size of tdTaudDOF incorrect" );

        // fill cells
        for ( uint iCol = 0; iCol < tNumSpaceDims; iCol++ )
        {
            aKijYjDOF( iCol + 1 ).set_size( tNumSpaceDims, tNumBases, 0.0 );
        }

        // clang-format off
            // put dof derivatives of stress tensor entries into right places
            if ( tNumSpaceDims == 2 )
            {
                // first column
                aKijYjDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 0, 0 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 2, 2 }, { 0, tNumBases - 1 } );

                // second column
                aKijYjDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 2, 2 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 2 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 1, 1 }, { 0, tNumBases - 1 } );
            }
            else if ( tNumSpaceDims == 3 )
            {
                // first column
                aKijYjDOF( 1 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 0, 0 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 1 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 5, 5 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 1 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 4, 4 }, { 0, tNumBases - 1 } );

                // second column
                aKijYjDOF( 2 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 5, 5 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 2 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 1, 1 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 2 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 3, 3 }, { 0, tNumBases - 1 } );

                // third column
                aKijYjDOF( 3 )( { 0, 0 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 4, 4 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 3 )( { 1, 1 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 3, 3 }, { 0, tNumBases - 1 } );
                aKijYjDOF( 3 )( { 2, 2 }, { 0, tNumBases - 1 } ) =
                        tdTaudDOF( { 2, 2 }, { 0, tNumBases - 1 } );
            }
            else
            {
                MORIS_ERROR( false,
                        "fn_FEM_IWG_Compressible_NS::eval_KijYjDOF_matrices - number of spatial dimensions != {2,3}" );
            }
        // clang-format on

        // temperature residual
        aKijYjDOF( tNumSpaceDims + 1 ) =
                aCM->dFluxdDOF( aDofType, CM_Function_Type::WORK ) - aCM->dFluxdDOF( aDofType, CM_Function_Type::THERMAL );
    }

    //------------------------------------------------------------------------------

    void eval_KijYjiDOF(
            const std::shared_ptr< Constitutive_Model > &aCM,
            Field_Interpolator_Manager                  *aLeaderFIManager,
            const Vector< Vector< MSI::Dof_Type > >     &aResidualDofTypes,
            const Vector< MSI::Dof_Type >               &aDofType,
            Matrix< DDRMat >                            &aKijYjiDOF )
    {

        // check inputs
        MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ),
                "fn_FEM_IWG_Compressible_NS::evaeval_KijYji_A_matrices - list of aResidualDofTypes not supported, see error messages above." );

        // get the velocity FI
        Field_Interpolator *tFIVelocity = aLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::VX );

        // get number of Space dimensions
        uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

        // get number of bases for the elements used
        uint tNumBases = tFIVelocity->get_number_of_space_time_bases();

        // initialize KijYj
        aKijYjiDOF.set_size( tNumSpaceDims + 2, tNumBases, 0.0 );

        // velocity residual
        aKijYjiDOF( { 1, tNumSpaceDims }, { 0, tNumBases - 1 } ) =
                aCM->ddivfluxdu( aDofType, CM_Function_Type::MECHANICAL ).matrix_data();

        // temperature residual
        aKijYjiDOF( { tNumSpaceDims + 1, tNumSpaceDims + 1 }, { 0, tNumBases - 1 } ) =
                aCM->ddivfluxdu( aDofType, CM_Function_Type::WORK ) - aCM->ddivfluxdu( aDofType, CM_Function_Type::THERMAL );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
