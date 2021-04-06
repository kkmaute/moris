/*
 * cl_FEM_IWG_Compressible_NS_Flux_Matrices_Var_Derivs.cpp
 *
 *  Created on: Mar 31, 2021
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

        void eval_VL_dAdY( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                const Matrix< DDRMat >                & aVL,
                moris::Cell< Matrix< DDRMat > >       & aVLdAdY )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - list of aResidualDofTypes not supported, see messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tP   = aMM->pressure()( 0 );
            real tT   = aMM->temperature()( 0 );
            real tRho = aMM->density()( 0 );
            real tR   = tP / tRho / tT;
            real tCv  = aMM->Cv()( 0 ); 
            real tU1  = tFIVelocity->val()( 0 );
            real tU2  = tFIVelocity->val()( 1 );

            // help values   
            real tU1sq = tU1*tU1;
            real tU2sq = tU2*tU2;
            real tC1 = 1.0/(tR*tT);
            real tC2 = tP/(tR*tT);
            real tC3 = 1.0/(tR*tT*tT);
            real tC4 = tP/(tR*tT*tT);

            // reset A matrices
            Matrix< DDRMat > tEmptyA( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
            aVLdAdY.assign( tNumSpaceDims + 2, tEmptyA );

            // check the pre-multiplication vector
            MORIS_ASSERT( aVL.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - length of pre-multiplication vector incorrect." );  

            // get values form pre-multiplication vector
            real tVL1 = aVL( 0 );
            real tVL2 = aVL( 1 );
            real tVL3 = aVL( 2 );
            real tVL4 = aVL( 3 );      

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // velocity magintude
                    //real tQ  = ( tU1sq + tU2sq ) / 2.0;

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tU3 = tFIVelocity->val()( 2 );
                    real tU3sq = tU3*tU3; 

                    // get last value from pre-multiplication vector
                    real tVL5 = aVL( 4 ); 

                    // velocity magintude
                    real tQ  = ( tU1sq + tU2sq + tU3sq ) / 2.0;

                    // evaluate A_0
                    aVLdAdY( 0 ) = { 
                        {                                                    0.0,  tC1*(tVL2 + tU1*tVL5),  tC1*(tVL3 + tU2*tVL5),  tC1*(tVL4 + tU3*tVL5),           -tC3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                        {                                  tC1*(tVL2 + tU1*tVL5),               tC2*tVL5,                    0.0,                    0.0,                                           -tC4*(tVL2 + tU1*tVL5) },
                        {                                  tC1*(tVL3 + tU2*tVL5),                    0.0,               tC2*tVL5,                    0.0,                                           -tC4*(tVL3 + tU2*tVL5) },
                        {                                  tC1*(tVL4 + tU3*tVL5),                    0.0,                    0.0,               tC2*tVL5,                                           -tC4*(tVL4 + tU3*tVL5) },
                        { -tC3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*(tVL2 + tU1*tVL5), -tC4*(tVL3 + tU2*tVL5), -tC4*(tVL4 + tU3*tVL5), (tC4*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                    // evaluate A_1
                    aVLdAdY( 1 ) = { 
                        {                                                                                             0.0, tC1*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),  tC1*tU1*(tVL3 + tU2*tVL5),  tC1*tU1*(tVL4 + tU3*tVL5),              -tC3*tU1*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                        { tC1*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),                                                                   tC2*(2.0*tVL2 + 3.0*tU1*tVL5),      tC2*(tVL3 + tU2*tVL5),      tC2*(tVL4 + tU3*tVL5), -tC4*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4) },
                        {                                                                       tC1*tU1*(tVL3 + tU2*tVL5),                                                                           tC2*(tVL3 + tU2*tVL5),               tC2*tU1*tVL5,                        0.0,                                              -tC4*tU1*(tVL3 + tU2*tVL5) },
                        {                                                                       tC1*tU1*(tVL4 + tU3*tVL5),                                                                           tC2*(tVL4 + tU3*tVL5),                        0.0,               tC2*tU1*tVL5,                                              -tC4*tU1*(tVL4 + tU3*tVL5) },
                        {                                      -tC3*tU1*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4),                         -tC4*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4), -tC4*tU1*(tVL3 + tU2*tVL5), -tC4*tU1*(tVL4 + tU3*tVL5),      (tC4*tU1*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                    // evaluate A_2
                    aVLdAdY( 2 ) = { 
                        {                                                                                             0.0,  tC1*tU2*(tVL2 + tU1*tVL5), tC1*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),  tC1*tU2*(tVL4 + tU3*tVL5),              -tC3*tU2*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                        {                                                                       tC1*tU2*(tVL2 + tU1*tVL5),               tC2*tU2*tVL5,                                                                           tC2*(tVL2 + tU1*tVL5),                        0.0,                                              -tC4*tU2*(tVL2 + tU1*tVL5) },
                        { tC1*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),      tC2*(tVL2 + tU1*tVL5),                                                                   tC2*(2.0*tVL3 + 3.0*tU2*tVL5),      tC2*(tVL4 + tU3*tVL5), -tC4*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4) },
                        {                                                                       tC1*tU2*(tVL4 + tU3*tVL5),                        0.0,                                                                           tC2*(tVL4 + tU3*tVL5),               tC2*tU2*tVL5,                                              -tC4*tU2*(tVL4 + tU3*tVL5) },
                        {                                      -tC3*tU2*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*tU2*(tVL2 + tU1*tVL5),                         -tC4*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4), -tC4*tU2*(tVL4 + tU3*tVL5),      (tC4*tU2*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                    // evaluate A_3
                    aVLdAdY( 3 ) = { 
                        {                                                                                             0.0,  tC1*tU3*(tVL2 + tU1*tVL5),  tC1*tU3*(tVL3 + tU2*tVL5), tC1*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tVL5/tC1 + tCv*tT*tVL5),              -tC3*tU3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                        {                                                                       tC1*tU3*(tVL2 + tU1*tVL5),               tC2*tU3*tVL5,                        0.0,                                                                           tC2*(tVL2 + tU1*tVL5),                                              -tC4*tU3*(tVL2 + tU1*tVL5) },
                        {                                                                       tC1*tU3*(tVL3 + tU2*tVL5),                        0.0,               tC2*tU3*tVL5,                                                                           tC2*(tVL3 + tU2*tVL5),                                              -tC4*tU3*(tVL3 + tU2*tVL5) },
                        { tC1*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tVL5/tC1 + tCv*tT*tVL5),      tC2*(tVL2 + tU1*tVL5),      tC2*(tVL3 + tU2*tVL5),                                                                   tC2*(2.0*tVL4 + 3.0*tU3*tVL5), -tC4*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3) },
                        {                                      -tC3*tU3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*tU3*(tVL2 + tU1*tVL5), -tC4*tU3*(tVL3 + tU2*tVL5),                         -tC4*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3),      (tC4*tU3*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_dAdY_VR( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                const Matrix< DDRMat >                & aVR,
                moris::Cell< Matrix< DDRMat > >       & adAdYVR )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::eval_dAdY_VR - list of aResidualDofTypes not supported, see messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            real tP   = aMM->pressure()( 0 );
            real tT   = aMM->temperature()( 0 );
            real tRho = aMM->density()( 0 );
            real tR   = tP / tRho / tT;
            real tCv  = aMM->Cv()( 0 ); 
            real tU1  = tFIVelocity->val()( 0 );
            real tU2  = tFIVelocity->val()( 1 );

            // help values   
            real tU1sq = tU1*tU1;
            real tU2sq = tU2*tU2;
            real tTsq = tT*tT;
            real tC3 = 1.0/(tR*tT*tT);

            // reset A matrices
            Matrix< DDRMat > tEmptyA( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
            adAdYVR.assign( tNumSpaceDims + 2, tEmptyA );

            // check the pre-multiplication vector
            MORIS_ASSERT( aVR.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_dAdY_VR - length of pre-multiplication vector incorrect." );  

            // get values form pre-multiplication vector
            real tVR1 = aVR( 0 );
            real tVR2 = aVR( 1 );
            real tVR3 = aVR( 2 );
            real tVR4 = aVR( 3 );      

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // velocity magintude
                    //real tQ  = ( tU1sq + tU2sq ) / 2.0;

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get the z-velocity
                    real tU3 = tFIVelocity->val()( 2 );
                    real tU3sq = tU3*tU3; 

                    // get last value from pre-multiplication vector
                    real tVR5 = aVR( 4 ); 

                    // velocity magintude
                    real tQ  = ( tU1sq + tU2sq + tU3sq ) / 2.0;

                    // evaluate A_0
                    adAdYVR( 0 ) = { 
                        {                      -tC3*tVR5,                             tC3*(tT*tVR2 - tU1*tVR5),                             tC3*(tT*tVR3 - tU2*tVR5),                             tC3*(tT*tVR4 - tU3*tVR5),                                  tC3*(tT*tU1*tVR2 - tQ*tVR5 + tT*tU2*tVR3 + tT*tU3*tVR4) },
                        {                              0,                             -tC3*(tP*tVR5 - tT*tVR1),                                                    0,                                                    0,                                             tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1) },
                        {                              0,                                                    0,                             -tC3*(tP*tVR5 - tT*tVR1),                                                    0,                                             tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1) },
                        {                              0,                                                    0,                                                    0,                             -tC3*(tP*tVR5 - tT*tVR1),                                             tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1) },
                        { (tC3*(2.0*tP*tVR5 - tT*tVR1))/tT, -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT, -(tC3*(tP*tT*tVR3 - 2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT, -(tC3*(tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT, -(tC3*(tQ*tT*tVR1 - 2.0*tP*tQ*tVR5 + tP*tT*tU1*tVR2 + tP*tT*tU2*tVR3 + tP*tT*tU3*tVR4))/tT } };

                    // evaluate A_1
                    adAdYVR( 1 ) = { 
                        {                             tC3*(tT*tVR2 - tU1*tVR5),                             tC3*tU1*(2.0*tT*tVR2 - tU1*tVR5),                                    tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR5),                                    tC3*(tT*tU1*tVR4 + tT*tU3*tVR2 - tU1*tU3*tVR5),                                   (tC3*((2.0*tVR2)/tC3 + 2.0*tCv*tTsq*tVR2 + 2.0*tQ*tT*tVR2 - 2.0*tQ*tU1*tVR5 + 2.0*tT*tU1sq*tVR2 + 2.0*tT*tU1*tU2*tVR3 + 2.0*tT*tU1*tU3*tVR4))/2.0 },
                        {                             -tC3*(tP*tVR5 - tT*tVR1),             2.0*tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                      tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                      tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1), (tC3*((2.0*tVR1)/tC3 + 2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR5 + 2.0*tQ*tT*tVR1 - 2.0*tP*tU1sq*tVR5 + 2.0*tT*tU1sq*tVR1 + 6.0*tP*tT*tU1*tVR2 + 2.0*tP*tT*tU2*tVR3 + 2.0*tP*tT*tU3*tVR4))/2.0 },
                        {                                                    0,                                                          0,                                      tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                                                                 0,                                                                                           tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1) },
                        {                                                    0,                                                          0,                                                                                 0,                                      tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                                                                           tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1) },
                        { -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT, -(tC3*tU1*(2.0*tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT, -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - 2.0*tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1))/tT, -(tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - 2.0*tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1))/tT,                                          -(tC3*(tP*tQ*tT*tVR2 - 2.0*tP*tQ*tU1*tVR5 + tP*tT*tU1sq*tVR2 + tQ*tT*tU1*tVR1 + tP*tT*tU1*tU2*tVR3 + tP*tT*tU1*tU3*tVR4))/tT } };
 
                    // evaluate A_2
                    adAdYVR( 2 ) = { 
                        {                             tC3*(tT*tVR3 - tU2*tVR5),                                    tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR5),                             tC3*tU2*(2.0*tT*tVR3 - tU2*tVR5),                                    tC3*(tT*tU2*tVR4 + tT*tU3*tVR3 - tU2*tU3*tVR5),                                   (tC3*((2.0*tVR3)/tC3 + 2.0*tCv*tTsq*tVR3 + 2.0*tQ*tT*tVR3 - 2.0*tQ*tU2*tVR5 + 2.0*tT*tU2sq*tVR3 + 2.0*tT*tU1*tU2*tVR2 + 2.0*tT*tU2*tU3*tVR4))/2.0 },
                        {                                                    0,                                      tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                                          0,                                                                                 0,                                                                                           tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1) },
                        {                             -tC3*(tP*tVR5 - tT*tVR1),                                      tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),             2.0*tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                      tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1), (tC3*((2.0*tVR1)/tC3 + 2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR5 + 2.0*tQ*tT*tVR1 - 2.0*tP*tU2sq*tVR5 + 2.0*tT*tU2sq*tVR1 + 2.0*tP*tT*tU1*tVR2 + 6.0*tP*tT*tU2*tVR3 + 2.0*tP*tT*tU3*tVR4))/2.0 },
                        {                                                    0,                                                                                 0,                                                          0,                                      tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                                                                           tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1) },
                        { -(tC3*(tP*tT*tVR3 - 2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT, -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - 2.0*tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1))/tT, -(tC3*tU2*(2.0*tP*tT*tVR3 - 2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT, -(tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - 2.0*tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1))/tT,                                          -(tC3*(tP*tQ*tT*tVR3 - 2.0*tP*tQ*tU2*tVR5 + tP*tT*tU2sq*tVR3 + tQ*tT*tU2*tVR1 + tP*tT*tU1*tU2*tVR2 + tP*tT*tU2*tU3*tVR4))/tT } };
 
                    // evaluate A_3
                    adAdYVR( 3 ) = { 
                        {                             tC3*(tT*tVR4 - tU3*tVR5),                                    tC3*(tT*tU1*tVR4 + tT*tU3*tVR2 - tU1*tU3*tVR5),                                    tC3*(tT*tU2*tVR4 + tT*tU3*tVR3 - tU2*tU3*tVR5),                             tC3*tU3*(2.0*tT*tVR4 - tU3*tVR5),                                   (tC3*((2.0*tVR4)/tC3 + 2.0*tCv*tTsq*tVR4 + 2.0*tQ*tT*tVR4 - 2.0*tQ*tU3*tVR5 + 2.0*tT*tU3sq*tVR4 + 2.0*tT*tU1*tU3*tVR2 + 2.0*tT*tU2*tU3*tVR3))/2.0 },
                        {                                                    0,                                      tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                                                 0,                                                          0,                                                                                           tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1) },
                        {                                                    0,                                                                                 0,                                      tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                          0,                                                                                           tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1) },
                        {                             -tC3*(tP*tVR5 - tT*tVR1),                                      tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                      tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),             2.0*tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1), (tC3*((2.0*tVR1)/tC3 + 2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR5 + 2.0*tQ*tT*tVR1 - 2.0*tP*tU3sq*tVR5 + 2.0*tT*tU3sq*tVR1 + 2.0*tP*tT*tU1*tVR2 + 2.0*tP*tT*tU2*tVR3 + 6.0*tP*tT*tU3*tVR4))/2.0 },
                        { -(tC3*(tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT, -(tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - 2.0*tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1))/tT, -(tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - 2.0*tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1))/tT, -(tC3*tU3*(2.0*tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT,                                          -(tC3*(tP*tQ*tT*tVR4 - 2.0*tP*tQ*tU3*tVR5 + tP*tT*tU3sq*tVR4 + tQ*tT*tU3*tVR1 + tP*tT*tU1*tU3*tVR2 + tP*tT*tU2*tU3*tVR3))/tT } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dAdY_VR() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_VL_dKdY( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                const Matrix< DDRMat >                         & aVL,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & aVLdKdY )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

            // initialize cells
            aVLdKdY.resize( tNumSpaceDims );
            for ( uint i = 0; i < tNumSpaceDims; i++ )
            {
                aVLdKdY( i ).resize( tNumSpaceDims );
            }

            // check the pre-multiplication vector
            MORIS_ASSERT( aVL.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - length of pre-multiplication vector incorrect." );      

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get last value from pre-multiplication vector
                    //real tVL4 = aVL( 3 ); 


                    break;
                }

                // for 3D
                case 3 :
                {
                    // get last value from pre-multiplication vector
                    real tVL5 = aVL( 4 ); 

                    // evaluate K11
                    aVLdKdY( 0 )( 0 ) = { 
                        { 0.0,                  0.0,      0.0,      0.0, 0.0 },
                        { 0.0, tVL5*(tLa + 2.0*tMu),      0.0,      0.0, 0.0 },
                        { 0.0,                  0.0, tMu*tVL5,      0.0, 0.0 },
                        { 0.0,                  0.0,      0.0, tMu*tVL5, 0.0 },
                        { 0.0,                  0.0,      0.0,      0.0, 0.0 } };

                    // evaluate K12
                    aVLdKdY( 0 )( 1 ) = { 
                        { 0.0, 0.0,      0.0, 0.0, 0.0 },
                        { 0.0, 0.0, tMu*tVL5, 0.0, 0.0 },
                        { 0.0, 0.0, tLa*tVL5, 0.0, 0.0 },
                        { 0.0, 0.0,      0.0, 0.0, 0.0 },
                        { 0.0, 0.0,      0.0, 0.0, 0.0 } };
 
                    // evaluate K13
                    aVLdKdY( 0 )( 2 ) = { 
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0, tMu*tVL5, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0, tLa*tVL5, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 } };

                    // ======= //

                    // evaluate K21
                    aVLdKdY( 1 )( 0 ) = { 
                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0, tLa*tVL5, 0.0, 0.0 },
                        { 0.0, tMu*tVL5,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0,      0.0, 0.0, 0.0 } };

                    // evaluate K22
                    aVLdKdY( 1 )( 1 ) = { 
                        { 0.0,      0.0,                  0.0,      0.0, 0.0 },
                        { 0.0, tMu*tVL5,                  0.0,      0.0, 0.0 },
                        { 0.0,      0.0, tVL5*(tLa + 2.0*tMu),      0.0, 0.0 },
                        { 0.0,      0.0,                  0.0, tMu*tVL5, 0.0 },
                        { 0.0,      0.0,                  0.0,      0.0, 0.0 } };
 
                    // evaluate K23
                    aVLdKdY( 1 )( 2 ) = { 
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0, tMu*tVL5, 0.0 },
                        { 0.0, 0.0, tLa*tVL5,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 } };

                    // ======= //
 
                     // evaluate K31
                    aVLdKdY( 2 )( 0 ) = { 
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0, tLa*tVL5, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0, tMu*tVL5, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 } };

                    // evaluate K32
                    aVLdKdY( 2 )( 1 ) = { 
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0, tLa*tVL5, 0.0 },
                        { 0.0, 0.0, tMu*tVL5,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 } };
 
                    // evaluate K33
                    aVLdKdY( 2 )( 2 ) = { 
                        { 0.0,      0.0,      0.0,                  0.0, 0.0 },
                        { 0.0, tMu*tVL5,      0.0,                  0.0, 0.0 },
                        { 0.0,      0.0, tMu*tVL5,                  0.0, 0.0 },
                        { 0.0,      0.0,      0.0, tVL5*(tLa + 2.0*tMu), 0.0 },
                        { 0.0,      0.0,      0.0,                  0.0, 0.0 } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_dKdY_VR( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                const Matrix< DDRMat >                         & aVR,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & adKdYVR )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

            // initialize cells
            adKdYVR.resize( tNumSpaceDims );
            for ( uint i = 0; i < tNumSpaceDims; i++ )
            {
                adKdYVR( i ).resize( tNumSpaceDims );
            }

            // check the pre-multiplication vector
            MORIS_ASSERT( aVR.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - length of pre-multiplication vector incorrect." );      

            // get values form pre-multiplication vector
            real tVR2 = aVR( 1 );
            real tVR3 = aVR( 2 );
            real tVR4 = aVR( 3 );    

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get last value from pre-multiplication vector
                    //real tVL4 = aVL( 3 ); 


                    break;
                }

                // for 3D
                case 3 :
                {
                    // evaluate K11
                    adKdYVR( 0 )( 0 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tVR2*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };

                    // evaluate K12
                    adKdYVR( 0 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR3 + tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 } };

                    // evaluate K13
                    adKdYVR( 0 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // ======= //

                    // evaluate K21
                    adKdYVR( 1 )( 0 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // evaluate K22
                    adKdYVR( 1 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0, tVR3*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };
 
                    // evaluate K23
                    adKdYVR( 1 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // ======= //
 
                     // evaluate K31
                    adKdYVR( 2 )( 0 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // evaluate K32
                    adKdYVR( 2 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };
 
                    // evaluate K33
                    adKdYVR( 2 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0, tVR4*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------
        //------------------------------------------------------------------------------

        void eval_VL_d2KdxdY( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                const Matrix< DDRMat >                         & aVL,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & aVLd2KdxdY)
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

            // initialize cells
            aVLd2KdxdY.resize( tNumSpaceDims );
            for ( uint i = 0; i < tNumSpaceDims; i++ )
            {
                aVLd2KdxdY( i ).resize( tNumSpaceDims + 1 );
            }

            // check the pre-multiplication vector
            MORIS_ASSERT( aVL.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_VL_d2KdxdY - length of pre-multiplication vector incorrect." );      

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get last value from pre-multiplication vector
                    //real tVL4 = aVL( 3 ); 

                    break;
                }

                // for 3D
                case 3 :
                {
                    // get last value from pre-multiplication vector
                    real tVL5 = aVL( 4 ); 

                    // evaluate Ki1,i - for Y,x
                    aVLd2KdxdY( 0 )( 1 ) = { 
                        { 0.0,                  0.0,      0.0,      0.0, 0.0 },
                        { 0.0, tVL5*(tLa + 2.0*tMu),      0.0,      0.0, 0.0 },
                        { 0.0,                  0.0, tMu*tVL5,      0.0, 0.0 },
                        { 0.0,                  0.0,      0.0, tMu*tVL5, 0.0 },
                        { 0.0,                  0.0,      0.0,      0.0, 0.0 } };

                    // evaluate Ki1,i - for Y,y
                    aVLd2KdxdY( 0 )( 2 ) = { 
                        { 0.0, 0.0,      0.0, 0.0, 0.0 },
                        { 0.0, 0.0, tMu*tVL5, 0.0, 0.0 },
                        { 0.0, 0.0, tLa*tVL5, 0.0, 0.0 },
                        { 0.0, 0.0,      0.0, 0.0, 0.0 },
                        { 0.0, 0.0,      0.0, 0.0, 0.0 } };
 
                    // evaluate Ki1,i - for Y,z
                    aVLd2KdxdY( 0 )( 3 ) = { 
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0, tMu*tVL5, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0, tLa*tVL5, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 } };

                    // ======= //

                    // evaluate Ki2,i - for Y,x
                    aVLd2KdxdY( 1 )( 1 ) = { 
                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0, tLa*tVL5, 0.0, 0.0 },
                        { 0.0, tMu*tVL5,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                        { 0.0,      0.0,      0.0, 0.0, 0.0 } };

                    // evaluate Ki2,i - for Y,y
                    aVLd2KdxdY( 1 )( 2 ) = { 
                        { 0.0,      0.0,                  0.0,      0.0, 0.0 },
                        { 0.0, tMu*tVL5,                  0.0,      0.0, 0.0 },
                        { 0.0,      0.0, tVL5*(tLa + 2.0*tMu),      0.0, 0.0 },
                        { 0.0,      0.0,                  0.0, tMu*tVL5, 0.0 },
                        { 0.0,      0.0,                  0.0,      0.0, 0.0 } };
 
                    // evaluate Ki2,i - for Y,z
                    aVLd2KdxdY( 1 )( 3 ) = { 
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0, tMu*tVL5, 0.0 },
                        { 0.0, 0.0, tLa*tVL5,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 } };

                    // ======= //
 
                     // evaluate Ki3,i - for Y,x
                    aVLd2KdxdY( 2 )( 1 ) = { 
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0, tLa*tVL5, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                        { 0.0, tMu*tVL5, 0.0,      0.0, 0.0 },
                        { 0.0,      0.0, 0.0,      0.0, 0.0 } };

                    // evaluate Ki3,i - for Y,y
                    aVLd2KdxdY( 2 )( 2 ) = { 
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0, tLa*tVL5, 0.0 },
                        { 0.0, 0.0, tMu*tVL5,      0.0, 0.0 },
                        { 0.0, 0.0,      0.0,      0.0, 0.0 } };
 
                    // evaluate Ki3,i - for Y,z
                    aVLd2KdxdY( 2 )( 3 ) = { 
                        { 0.0,      0.0,      0.0,                  0.0, 0.0 },
                        { 0.0, tMu*tVL5,      0.0,                  0.0, 0.0 },
                        { 0.0,      0.0, tMu*tVL5,                  0.0, 0.0 },
                        { 0.0,      0.0,      0.0, tVL5*(tLa + 2.0*tMu), 0.0 },
                        { 0.0,      0.0,      0.0,                  0.0, 0.0 } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_d2KdxdY() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_d2KdxdY_VR( 
                std::shared_ptr< Property >                      aPropDynamicViscosity,  
                std::shared_ptr< Property >                      aPropThermalConductivity,
                Field_Interpolator_Manager                     * aMasterFIManager,
                const Matrix< DDRMat >                         & aVR,
                moris::Cell< moris::Cell< Matrix< DDRMat > > > & ad2KdxdYVR)
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

            // initialize cells
            ad2KdxdYVR.resize( tNumSpaceDims );
            for ( uint i = 0; i < tNumSpaceDims; i++ )
            {
                ad2KdxdYVR( i ).resize( tNumSpaceDims + 1 );
            }

            // check the pre-multiplication vector
            MORIS_ASSERT( aVR.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_d2KdxdY_VR - length of pre-multiplication vector incorrect." );      

            // get values form pre-multiplication vector
            real tVR2 = aVR( 1 );
            real tVR3 = aVR( 2 );
            real tVR4 = aVR( 3 );   

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get last value from pre-multiplication vector
                    //real tVL4 = aVL( 3 ); 

                    break;
                }

                // for 3D
                case 3 :
                {
                    // evaluate Ki1,i - for Y,x
                    ad2KdxdYVR( 0 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tVR2*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };

                    // evaluate Ki1,i - for Y,y
                    ad2KdxdYVR( 0 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR3 + tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 },
                        { 0.0, 0.0, 0.0, 0.0,                 0.0 } };
 
                    // evaluate Ki1,i - for Y,z
                    ad2KdxdYVR( 0 )( 3 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // ======= //

                    // evaluate Ki2,i - for Y,x
                    ad2KdxdYVR( 1 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // evaluate Ki2,i - for Y,y
                    ad2KdxdYVR( 1 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0, tVR3*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };
 
                    // evaluate Ki2,i - for Y,z
                    ad2KdxdYVR( 1 )( 3 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // ======= //
 
                     // evaluate Ki3,i - for Y,x
                    ad2KdxdYVR( 2 )( 1 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };

                    // evaluate Ki3,i - for Y,y
                    ad2KdxdYVR( 2 )( 2 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 },
                        { 0.0, 0.0, 0.0, 0.0, tMu*tVR4 },
                        { 0.0, 0.0, 0.0, 0.0, tLa*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0,      0.0 } };
 
                    // evaluate Ki3,i - for Y,z
                    ad2KdxdYVR( 2 )( 3 ) = { 
                        { 0.0, 0.0, 0.0, 0.0,                0.0 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR2 },
                        { 0.0, 0.0, 0.0, 0.0,           tMu*tVR3 },
                        { 0.0, 0.0, 0.0, 0.0, tVR4*(tLa+2.0*tMu) },
                        { 0.0, 0.0, 0.0, 0.0,                0.0 } };

                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_d2KdxdY_VR() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
