/*
 * fn_FEM_IWG_Compressible_NS_Flux_Matrices.cpp
 *
 *  Created on: Mar 18, 2021
 *      Author: wunsch
 */

#include "fn_FEM_IWG_Compressible_NS.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris
{
    namespace fem
    {
        //------------------------------------------------------------------------------

        void eval_A( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                moris::Cell< Matrix< DDRMat > >       & aAMats )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::eval_A_matrices - list of aResidualDofTypes not supported, see messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tRho = aMM->density()( 0 );
            real tUx = tFIVelocity->val()( 0 );
            real tUy = tFIVelocity->val()( 1 );
            real tAlphaP = aMM->AlphaP()( 0 );
            real tBetaT = aMM->BetaT()( 0 );
            real tEtot = aCM->Energy()( 0 );

            // constants computed as in M. Polner's 2005 PhD thesis
            real tE1p = tBetaT * tEtot;
            real tE2p = tE1p + 1.0;
            real tE3p = tEtot + aMM->pressure()( 0 );
            real tE4p = -1.0 * tAlphaP * tEtot + tRho * aMM->Cv()( 0 );         

            // reset A matrices
            Matrix< DDRMat > tEmptyA( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
            aAMats.assign( tNumSpaceDims + 2, tEmptyA );
            //mA.resize( tNumSpaceDims + 2 );          

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // evaluate A_0
                    aAMats( 0 ) = { 
                        {       tRho * tBetaT,        0.0,        0.0,       -tRho * tAlphaP }, 
                        { tUx * tRho * tBetaT,       tRho,        0.0, tUx * -tRho * tAlphaP }, 
                        { tUy * tRho * tBetaT,        0.0,       tRho, tUy * -tRho * tAlphaP }, 
                        {                tE1p, tUx * tRho, tUy * tRho,                  tE4p } };

                    // evaluate A_1
                    aAMats( 1 ) = { 
                        {             tRho * tBetaT * tUx,                    tRho,              0.0,       -tRho * tAlphaP * tUx }, 
                        { 1.0 + tUx * tRho * tBetaT * tUx,        2.0 * tUx * tRho,              0.0, tUx * -tRho * tAlphaP * tUx }, 
                        {       tUy * tRho * tBetaT * tUx,              tUy * tRho,       tUx * tRho, tUy * -tRho * tAlphaP * tUx },  
                        {                      tE2p * tUx, tE3p + tUx * tUx * tRho, tUy * tUx * tRho,                  tE4p * tUx } };

                    // evaluate A_2
                    aAMats( 2 ) = { 
                        {             tRho * tBetaT * tUy,              0.0,                    tRho,       -tRho * tAlphaP * tUy }, 
                        {       tUx * tRho * tBetaT * tUy,       tUy * tRho,              tUx * tRho, tUx * -tRho * tAlphaP * tUy }, 
                        { 1.0 + tUy * tRho * tBetaT * tUy,              0.0,        2.0 * tUy * tRho, tUy * -tRho * tAlphaP * tUy }, 
                        {                      tE2p * tUy, tUx * tUy * tRho, tE3p + tUy * tUy * tRho,                  tE4p * tUy } };

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tUz = tFIVelocity->val()( 2 );

                    // evaluate A_0
                    aAMats( 0 ) = { 
                        {       tRho * tBetaT,        0.0,        0.0,        0.0,       -tRho * tAlphaP }, 
                        { tUx * tRho * tBetaT,       tRho,        0.0,        0.0, tUx * -tRho * tAlphaP }, 
                        { tUy * tRho * tBetaT,        0.0,       tRho,        0.0, tUy * -tRho * tAlphaP }, 
                        { tUz * tRho * tBetaT,        0.0,        0.0,       tRho, tUz * -tRho * tAlphaP }, 
                        {                tE1p, tUx * tRho, tUy * tRho, tUz * tRho,                  tE4p } };      

                    // evaluate A_1
                    aAMats( 1 ) = { 
                        {             tRho * tBetaT * tUx,                    tRho,              0.0,              0.0,       -tRho * tAlphaP * tUx }, 
                        { 1.0 + tUx * tRho * tBetaT * tUx,        2.0 * tUx * tRho,              0.0,              0.0, tUx * -tRho * tAlphaP * tUx }, 
                        {       tUy * tRho * tBetaT * tUx,              tUy * tRho,       tUx * tRho,              0.0, tUy * -tRho * tAlphaP * tUx }, 
                        {       tUz * tRho * tBetaT * tUx,              tUz * tRho,              0.0,       tUx * tRho, tUz * -tRho * tAlphaP * tUx }, 
                        {                      tE2p * tUx, tE3p + tUx * tUx * tRho, tUy * tUx * tRho, tUz * tUx * tRho,                  tE4p * tUx } };

                    // evaluate A_2
                    aAMats( 2 ) = { 
                        {             tRho * tBetaT * tUy,              0.0,                    tRho,              0.0,       -tRho * tAlphaP * tUy }, 
                        {       tUx * tRho * tBetaT * tUy,       tUy * tRho,              tUx * tRho,              0.0, tUx * -tRho * tAlphaP * tUy }, 
                        { 1.0 + tUy * tRho * tBetaT * tUy,              0.0,        2.0 * tUy * tRho,              0.0, tUy * -tRho * tAlphaP * tUy }, 
                        {       tUz * tRho * tBetaT * tUy,              0.0,              tUz * tRho,       tUy * tRho, tUz * -tRho * tAlphaP * tUy }, 
                        {                      tE2p * tUy, tUx * tUy * tRho, tE3p + tUy * tUy * tRho, tUz * tUy * tRho,                  tE4p * tUy } };


                    // evaluate A_3
                    aAMats( 3 ) = { 
                        {             tRho * tBetaT * tUz,              0.0,              0.0,                    tRho,       -tRho * tAlphaP * tUz }, 
                        {       tUx * tRho * tBetaT * tUz,       tUz * tRho,              0.0,              tUx * tRho, tUx * -tRho * tAlphaP * tUz }, 
                        {       tUy * tRho * tBetaT * tUz,              0.0,       tUz * tRho,              tUy * tRho, tUy * -tRho * tAlphaP * tUz }, 
                        { 1.0 + tUz * tRho * tBetaT * tUz,              0.0,              0.0,        2.0 * tUz * tRho, tUz * -tRho * tAlphaP * tUz }, 
                        {                      tE2p * tUz, tUx * tUz * tRho, tUy * tUz * tRho, tE3p + tUz * tUz * tRho,                  tE4p * tUz } };                        

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::assemble_A_matrices() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_K( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & aKMats )
        {  
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tKa = aPropThermalConductivity->val()( 0 );
            real tMu = aPropDynamicViscosity->val()( 0 );
            real tLa = -2.0 * tMu / 3.0;
            real tCh = 2.0 * tMu + tLa;
            real tUx = tFIVelocity->val()( 0 );
            real tUy = tFIVelocity->val()( 1 );

            // set number of K matrices
            aKMats.resize( tNumSpaceDims );
            for ( uint iDim = 0; iDim < tNumSpaceDims ; iDim++)
            {
                aKMats( iDim ).resize( tNumSpaceDims );
            }            

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // evaluate K_11
                    aKMats( 0 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tCh,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tCh, tUy * tMu, tKa } };

                    // evaluate K_12
                    aKMats( 0 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tLa, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0, tUy * tMu, tUx * tLa, 0.0 } };

                    // =======================

                    // evaluate K_21
                    aKMats( 1 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0,       tLa,       0.0, 0.0 }, 
                        { 0.0, tUy * tLa, tUx * tMu, 0.0 } };

                        // evaluate K_22
                    aKMats( 1 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tCh, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tCh, tKa } };

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tUz = tFIVelocity->val()( 2 );

                    // evaluate K_11
                    aKMats( 0 )( 0 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tCh,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tCh, tUy * tMu, tUz * tMu, tKa } };

                    // evaluate K_12
                    aKMats( 0 )( 1 ) = { 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       tLa, 0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0, tUy * tMu, tUx * tLa, 0.0, 0.0 } };

                    // evaluate K_13
                    aKMats( 0 )( 2 ) = { 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0, 0.0,       tLa, 0.0 }, 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu, 0.0,       0.0, 0.0 }, 
                        { 0.0, tUz * tMu, 0.0, tUx * tLa, 0.0 } };

                    // =======================

                    // evaluate K_21
                    aKMats( 1 )( 0 ) = { 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu, 0.0, 0.0 }, 
                        { 0.0,       tLa,       0.0, 0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0, 0.0, 0.0 }, 
                        { 0.0, tUy * tLa, tUx * tMu, 0.0, 0.0 } };

                        // evaluate K_22
                    aKMats( 1 )( 1 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tCh,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tCh, tUz * tMu, tKa } };

                        // evaluate K_23
                    aKMats( 1 )( 2 ) = { 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       tLa, 0.0 }, 
                        { 0.0, 0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0, 0.0, tUy * tMu, tUy * tLa, 0.0 } };

                    // =======================

                    // evaluate K_31
                    aKMats( 2 )( 0 ) = { 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0, 0.0,       tMu, 0.0 }, 
                        { 0.0,       0.0, 0.0,       0.0, 0.0 }, 
                        { 0.0,       tLa, 0.0,       0.0, 0.0 }, 
                        { 0.0, tUz * tLa, 0.0, tUx * tMu, 0.0 } };

                    // evaluate K_32
                    aKMats( 2 )( 1 ) = { 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0, 0.0,       0.0,       tMu, 0.0 }, 
                        { 0.0, 0.0,       tLa,       0.0, 0.0 }, 
                        { 0.0, 0.0, tUz * tLa, tUy * tMu, 0.0 } };

                    // evaluate K_33
                    aKMats( 2 )( 2 ) = { 
                        { 0.0,       0.0,       0.0,       0.0, 0.0 }, 
                        { 0.0,       tMu,       0.0,       0.0, 0.0 }, 
                        { 0.0,       0.0,       tMu,       0.0, 0.0 }, 
                        { 0.0,       0.0,       0.0,       tCh, 0.0 }, 
                        { 0.0, tUx * tMu, tUy * tMu, tUz * tCh, tKa } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::assemble_A_matrices() - Number of space dimensions must be 2 or 3" );
                };
            }
        }                    

        //------------------------------------------------------------------------------

        void eval_KijYj( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes,
                Matrix< DDRMat >                      & aKijYj )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::evaeval_KijYjl_A_matrices - list of aResidualDofTypes not supported, see error messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // initialize KijYj
            aKijYj.set_size( tNumSpaceDims, tNumSpaceDims + 2, 0.0 );

            // velocity term
            aKijYj( { 0, tNumSpaceDims - 1 }, { 1, tNumSpaceDims } ) = 
                    unfold_flat_tensor( aCM->flux( CM_Function_Type::WORK ) ).matrix_data();

            // temperature term
            aKijYj( { 0, tNumSpaceDims - 1 }, { tNumSpaceDims + 1, tNumSpaceDims + 1 } ) = 
                    aCM->flux( CM_Function_Type::WORK ) - aCM->flux( CM_Function_Type::THERMAL );
        }

        //------------------------------------------------------------------------------

        void eval_KijYji( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes,
                Matrix< DDRMat >                      & aKijYji )
        {

            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::evaeval_KijYji_A_matrices - list of aResidualDofTypes not supported, see error messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // initialize KijYj
            aKijYji.set_size( tNumSpaceDims + 2, 1, 0.0 );

            // velocity residual
            aKijYji( { 1, tNumSpaceDims } ) = 
                    aCM->divflux( CM_Function_Type::MECHANICAL ).matrix_data();

            // temperature residual
            aKijYji( { tNumSpaceDims + 1, tNumSpaceDims + 1 } ) =
                    aCM->divflux( CM_Function_Type::WORK ) - aCM->divflux( CM_Function_Type::THERMAL );
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void eval_KijYjDOF( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes,
                const moris::Cell< MSI::Dof_Type >    & aDofType,
                moris::Cell< Matrix< DDRMat > >       & aKijYjDOF )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::eval_KijYjDOF_matrices - list of aResidualDofTypes not supported, see error messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get number of bases for the elements used
            uint tNumBases = tFIVelocity->get_number_of_space_time_bases();  

            // initialize dKijYj/dDOF, each cell entry represents the dof derivs of one column of Kij*Y,j
            aKijYjDOF.reserve( tNumSpaceDims + 2 );
            
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

            // temperature residual
            aKijYjDOF( tNumSpaceDims + 1  ) = 
                    aCM->dFluxdDOF( aDofType, CM_Function_Type::WORK ) - aCM->dFluxdDOF( aDofType, CM_Function_Type::THERMAL );
        }

        //------------------------------------------------------------------------------

        void eval_KijYjiDOF( 
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes,
                const moris::Cell< MSI::Dof_Type >    & aDofType,
                Matrix< DDRMat >                      & aKijYjiDOF )
        {

            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::evaeval_KijYji_A_matrices - list of aResidualDofTypes not supported, see error messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
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

        Matrix< DDRMat > unfold_flat_tensor( const Matrix< DDRMat > & aFlattenedTensor )
        {
            // get length of the flattened tensor
            uint tLength = aFlattenedTensor.length();

            switch ( tLength )
            {
                // 2D: convert 3x1 flattened tensor back to 2x2 matrix
                case 3 :
                {
                    return {
                            { aFlattenedTensor( 0 ), aFlattenedTensor( 2 ) },
                            { aFlattenedTensor( 2 ), aFlattenedTensor( 1 ) } };
                    break;
                }
 
                // 3D: convert 6x1 flattened tensor back to 3x3 matrix
                case 6 :
                {
                    return {
                            { aFlattenedTensor( 0 ), aFlattenedTensor( 5 ), aFlattenedTensor( 4 ) },
                            { aFlattenedTensor( 5 ), aFlattenedTensor( 1 ), aFlattenedTensor( 3 ) },
                            { aFlattenedTensor( 4 ), aFlattenedTensor( 3 ), aFlattenedTensor( 2 ) } };
                    break;
                }
                
                default:
                {
                    MORIS_ERROR( false, 
                            "fn_FEM_IWG_Compressible_NS::unfold_flat_tensor - aFlattenedTensor must be a vector of length 3 (2D) or 6 (3D)." );
                    return { { 0.0 } };
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
