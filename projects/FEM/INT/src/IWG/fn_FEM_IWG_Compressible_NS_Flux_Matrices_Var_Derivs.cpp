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
                const uint                              aI,
                Matrix< DDRMat >                      & aVLdAdY )
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
                    real tQ  = ( tU1sq + tU2sq ) / 2.0;

                    // compute pre-multiplied A-matrix deriv for requested index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            aVLdAdY = {
                                {                                         0.0,  tC1*(tVL2 + tU1*tVL4),  tC1*(tVL3 + tU2*tVL4),                     -tC3*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3) },
                                {                       tC1*(tVL2 + tU1*tVL4),               tC2*tVL4,                    0.0,                                          -tC4*(tVL2 + tU1*tVL4) },
                                {                       tC1*(tVL3 + tU2*tVL4),                    0.0,               tC2*tVL4,                                          -tC4*(tVL3 + tU2*tVL4) },
                                { -tC3*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3), -tC4*(tVL2 + tU1*tVL4), -tC4*(tVL3 + tU2*tVL4), (tC4*(2.0*tVL1 + 2.0*tQ*tVL4 + 2.0*tU1*tVL2 + 2.0*tU2*tVL3))/tT } };
 
                            break; 
                        }

                        case 1 :
                        {
                            aVLdAdY = {
                                {                                                                                  0.0, tC1*(tVL4*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL4 + tU2*tVL3 + tVL4/tC1 + tCv*tT*tVL4),  tC1*tU1*(tVL3 + tU2*tVL4),                     -tC3*tU1*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3) },
                                { tC1*(tVL4*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL4 + tU2*tVL3 + tVL4/tC1 + tCv*tT*tVL4),                                                        tC2*(2.0*tVL2 + 3.0*tU1*tVL4),      tC2*(tVL3 + tU2*tVL4),        -tC4*(tVL4*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL4 + tU2*tVL3) },
                                {                                                            tC1*tU1*(tVL3 + tU2*tVL4),                                                                tC2*(tVL3 + tU2*tVL4),               tC2*tU1*tVL4,                                          -tC4*tU1*(tVL3 + tU2*tVL4) },
                                {                                      -tC3*tU1*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3),                         -tC4*(tVL4*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL4 + tU2*tVL3), -tC4*tU1*(tVL3 + tU2*tVL4), (tC4*tU1*(2.0*tVL1 + 2.0*tQ*tVL4 + 2.0*tU1*tVL2 + 2.0*tU2*tVL3))/tT } };
 
                            break; 
                        }

                        case 2 :
                        {
                            aVLdAdY = {
                                {                                                                                  0.0,  tC1*tU2*(tVL2 + tU1*tVL4), tC1*(tVL4*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL4 + tU1*tVL2 + tVL4/tC1 + tCv*tT*tVL4),                     -tC3*tU2*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3) },
                                {                                                            tC1*tU2*(tVL2 + tU1*tVL4),               tC2*tU2*tVL4,                                                                tC2*(tVL2 + tU1*tVL4),                                          -tC4*tU2*(tVL2 + tU1*tVL4) },
                                { tC1*(tVL4*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL4 + tU1*tVL2 + tVL4/tC1 + tCv*tT*tVL4),      tC2*(tVL2 + tU1*tVL4),                                                        tC2*(2.0*tVL3 + 3.0*tU2*tVL4),        -tC4*(tVL4*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL4 + tU1*tVL2) },
                                {                                      -tC3*tU2*(tVL1 + tQ*tVL4 + tU1*tVL2 + tU2*tVL3), -tC4*tU2*(tVL2 + tU1*tVL4),                         -tC4*(tVL4*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL4 + tU1*tVL2), (tC4*tU2*(2.0*tVL1 + 2.0*tQ*tVL4 + 2.0*tU1*tVL2 + 2.0*tU2*tVL3))/tT } };
 
                            break; 
                        }
                    
                        default:
                        {
                            // error, A-matrix index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - index for A-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices

                    // break for 2D case
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

                    // compute pre-multiplied A-matrix deriv for requested index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            aVLdAdY = {
                                {                                                    0.0,  tC1*(tVL2 + tU1*tVL5),  tC1*(tVL3 + tU2*tVL5),  tC1*(tVL4 + tU3*tVL5),           -tC3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                                {                                  tC1*(tVL2 + tU1*tVL5),               tC2*tVL5,                    0.0,                    0.0,                                           -tC4*(tVL2 + tU1*tVL5) },
                                {                                  tC1*(tVL3 + tU2*tVL5),                    0.0,               tC2*tVL5,                    0.0,                                           -tC4*(tVL3 + tU2*tVL5) },
                                {                                  tC1*(tVL4 + tU3*tVL5),                    0.0,                    0.0,               tC2*tVL5,                                           -tC4*(tVL4 + tU3*tVL5) },
                                { -tC3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*(tVL2 + tU1*tVL5), -tC4*(tVL3 + tU2*tVL5), -tC4*(tVL4 + tU3*tVL5), (tC4*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                            break; 
                        }

                        case 1 :
                        {
                            aVLdAdY = {
                                {                                                                                             0.0, tC1*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),  tC1*tU1*(tVL3 + tU2*tVL5),  tC1*tU1*(tVL4 + tU3*tVL5),              -tC3*tU1*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                                { tC1*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),                                                                   tC2*(2.0*tVL2 + 3.0*tU1*tVL5),      tC2*(tVL3 + tU2*tVL5),      tC2*(tVL4 + tU3*tVL5), -tC4*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4) },
                                {                                                                       tC1*tU1*(tVL3 + tU2*tVL5),                                                                           tC2*(tVL3 + tU2*tVL5),               tC2*tU1*tVL5,                        0.0,                                              -tC4*tU1*(tVL3 + tU2*tVL5) },
                                {                                                                       tC1*tU1*(tVL4 + tU3*tVL5),                                                                           tC2*(tVL4 + tU3*tVL5),                        0.0,               tC2*tU1*tVL5,                                              -tC4*tU1*(tVL4 + tU3*tVL5) },
                                {                                      -tC3*tU1*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4),                         -tC4*(tVL5*tU1sq + 2.0*tVL2*tU1 + tVL1 + tQ*tVL5 + tU2*tVL3 + tU3*tVL4), -tC4*tU1*(tVL3 + tU2*tVL5), -tC4*tU1*(tVL4 + tU3*tVL5),      (tC4*tU1*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                            break; 
                        }

                        case 2 :
                        {
                            aVLdAdY = {
                                {                                                                                             0.0,  tC1*tU2*(tVL2 + tU1*tVL5), tC1*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),  tC1*tU2*(tVL4 + tU3*tVL5),              -tC3*tU2*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                                {                                                                       tC1*tU2*(tVL2 + tU1*tVL5),               tC2*tU2*tVL5,                                                                           tC2*(tVL2 + tU1*tVL5),                        0.0,                                              -tC4*tU2*(tVL2 + tU1*tVL5) },
                                { tC1*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4 + tVL5/tC1 + tCv*tT*tVL5),      tC2*(tVL2 + tU1*tVL5),                                                                   tC2*(2.0*tVL3 + 3.0*tU2*tVL5),      tC2*(tVL4 + tU3*tVL5), -tC4*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4) },
                                {                                                                       tC1*tU2*(tVL4 + tU3*tVL5),                        0.0,                                                                           tC2*(tVL4 + tU3*tVL5),               tC2*tU2*tVL5,                                              -tC4*tU2*(tVL4 + tU3*tVL5) },
                                {                                      -tC3*tU2*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*tU2*(tVL2 + tU1*tVL5),                         -tC4*(tVL5*tU2sq + 2.0*tVL3*tU2 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU3*tVL4), -tC4*tU2*(tVL4 + tU3*tVL5),      (tC4*tU2*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                            break; 
                        }

                        case 3 :
                        {
                            aVLdAdY = { 
                                {                                                                                             0.0,  tC1*tU3*(tVL2 + tU1*tVL5),  tC1*tU3*(tVL3 + tU2*tVL5), tC1*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tVL5/tC1 + tCv*tT*tVL5),              -tC3*tU3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4) },
                                {                                                                       tC1*tU3*(tVL2 + tU1*tVL5),               tC2*tU3*tVL5,                        0.0,                                                                           tC2*(tVL2 + tU1*tVL5),                                              -tC4*tU3*(tVL2 + tU1*tVL5) },
                                {                                                                       tC1*tU3*(tVL3 + tU2*tVL5),                        0.0,               tC2*tU3*tVL5,                                                                           tC2*(tVL3 + tU2*tVL5),                                              -tC4*tU3*(tVL3 + tU2*tVL5) },
                                { tC1*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tVL5/tC1 + tCv*tT*tVL5),      tC2*(tVL2 + tU1*tVL5),      tC2*(tVL3 + tU2*tVL5),                                                                   tC2*(2.0*tVL4 + 3.0*tU3*tVL5), -tC4*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3) },
                                {                                      -tC3*tU3*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4), -tC4*tU3*(tVL2 + tU1*tVL5), -tC4*tU3*(tVL3 + tU2*tVL5),                         -tC4*(tVL5*tU3sq + 2.0*tVL4*tU3 + tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3),      (tC4*tU3*2.0*(tVL1 + tQ*tVL5 + tU1*tVL2 + tU2*tVL3 + tU3*tVL4))/tT } };

                            break;
                        }
                    
                        default:
                        {
                            // error, A-matrix index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - index for A-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices

                    // break for 3D case
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

        void eval_dAdY( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                const uint                              aYind,
                moris::Cell< Matrix< DDRMat > >       & adAdY )
        {
            // check inputs
            MORIS_ASSERT( check_residual_dof_types( aResidualDofTypes ), 
                    "fn_FEM_IWG_Compressible_NS::eval_dAdY - list of aResidualDofTypes not supported, see messages above." );

            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // size dAdY correctly
            adAdY.resize( tNumSpaceDims + 1 );

            // // get commonly used values
            real tP   = aMM->pressure()( 0 );
            real tT   = aMM->temperature()( 0 );
            real tRho = aMM->density()( 0 );
            real tR   = tP / tRho / tT;
            //real tCv  = aMM->Cv()( 0 ); 
            real tCp  = aMM->Cp()( 0 ); 
            real tEint= aMM->Eint()( 0 );
            real tU1  = tFIVelocity->val()( 0 );
            real tU2  = tFIVelocity->val()( 1 );

            // help values   
            real tU1sq = tU1*tU1;
            real tU2sq = tU2*tU2;
            real tC1 = 1.0/(tR*tT);
            real tC2 = tP/(tR*tT);
            real tC3 = 1.0/(tR*tT*tT);
            real tC4 = tP/(tR*tT*tT);

            // assemble matrices based on number of spatial dimensions
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get the 1/2 * v^2 term
                    real tQ = 0.5 * ( tU1sq * tU1sq + tU2sq * tU2sq );

                    // assemble matrices for requested state variable derivative
                    switch ( aYind )
                    {
                        // for PRESSURE
                        case 0 :
                        {
                            adAdY( 0 ) = {
                            	{ 0.0,     0.0,     0.0,                               -tC3 },
                            	{ 0.0,     tC1,     0.0,                           -tC3*tU1 },
                            	{ 0.0,     0.0,     tC1,                           -tC3*tU2 },
                            	{ 0.0, tC1*tU1, tC1*tU2, tC1*tCp - tC3*(tEint + tQ+1.0/tC1) } };
                             
                            adAdY( 1 ) = {
                            	{ 0.0,                                tC1,         0.0,                                 -tC3*tU1 },
                            	{ 0.0,                        2.0*tC1*tU1,         0.0,                               -tC3*tU1sq },
                            	{ 0.0,                            tC1*tU2,     tC1*tU1,                             -tC3*tU1*tU2 },
                            	{ 0.0, tC1*tEint + tC1*tQ + tC1*tU1sq+1.0, tC1*tU1*tU2, tU1*(tC1*tCp - tC3*(tEint + tQ+1.0/tC1)) } };                             
                             
                            adAdY( 2 ) = {
                            	{ 0.0,         0.0,                                tC1,                                 -tC3*tU2 },
                            	{ 0.0,     tC1*tU2,                            tC1*tU1,                             -tC3*tU1*tU2 },
                            	{ 0.0,         0.0,                        2.0*tC1*tU2,                               -tC3*tU2sq },
                            	{ 0.0, tC1*tU1*tU2, tC1*tEint + tC1*tQ + tC1*tU2sq+1.0, tU2*(tC1*tCp - tC3*(tEint + tQ+1.0/tC1)) } };

                            // break for PRESSURE case
                            break;
                        }

                        // for VX
                        case 1 :
                        {
                            adAdY( 0 ) = {
                            	{     0.0, 0.0, 0.0,        0 },
                            	{     tC1, 0.0, 0.0,     -tC4 },
                            	{     0.0, 0.0, 0.0,        0 },
                            	{ tC1*tU1, tC2, 0.0, -tC4*tU1 } };                            
                             
                            adAdY( 1 ) = {
                            	{                                tC1,         0.0,     0.0,                                           -tC4 },
                            	{                        2.0*tC1*tU1,     2.0*tC2,     0.0,                                   -2.0*tC4*tU1 },
                            	{                            tC1*tU2,         0.0,     tC2,                                       -tC4*tU2 },
                            	{ tC1*tEint + tC1*tQ + tC1*tU1sq+1.0, 3.0*tC2*tU1, tC2*tU2, tC2*tCp - tC4*tU1sq - tC4*(tEint + tQ+1.0/tC1) } };                            
                             
                            adAdY( 2 ) = {
                            	{         0.0,     0.0,     0.0,            0 },
                            	{     tC1*tU2,     0.0,     tC2,     -tC4*tU2 },
                            	{         0.0,     0.0,     0.0,            0 },
                            	{ tC1*tU1*tU2, tC2*tU2, tC2*tU1, -tC4*tU1*tU2 } };
                            
                            // break for VX case
                            break;
                        }

                        // for VY
                        case 2 :
                        {
                            adAdY( 0 ) = {
                            	{     0.0, 0.0, 0.0,        0 },
                            	{     0.0, 0.0, 0.0,        0 },
                            	{     tC1, 0.0, 0.0,     -tC4 },
                            	{ tC1*tU2, 0.0, tC2, -tC4*tU2 } };                           
                             
                            adAdY( 1 ) = {
                            	{         0.0,     0.0,     0.0,            0 },
                            	{         0.0,     0.0,     0.0,            0 },
                            	{     tC1*tU1,     tC2,     0.0,     -tC4*tU1 },
                            	{ tC1*tU1*tU2, tC2*tU2, tC2*tU1, -tC4*tU1*tU2 } };                        
                             
                            adAdY( 2 ) = {
                            	{                                tC1,     0.0,         0.0,                                           -tC4 },
                            	{                            tC1*tU1,     tC2,         0.0,                                       -tC4*tU1 },
                            	{                        2.0*tC1*tU2,     0.0,     2.0*tC2,                                   -2.0*tC4*tU2 },
                            	{ tC1*tEint + tC1*tQ + tC1*tU2sq+1.0, tC2*tU1, 3.0*tC2*tU2, tC2*tCp - tC4*tU2sq - tC4*(tEint + tQ+1.0/tC1) } };
                            
                            // break for VY case
                            break;
                        }

                        // for TEMP
                        case 3 :
                        {
                            adAdY( 0 ) = {
                            	{                               -tC3,      0.0,      0.0,                                  (2.0*tC4)/tT },
                            	{                           -tC3*tU1,     -tC4,      0.0,                              (2.0*tC4*tU1)/tT },
                            	{                           -tC3*tU2,      0.0,     -tC4,                              (2.0*tC4*tU2)/tT },
                            	{ tC1*tCp - tC3*(tEint + tQ+1.0/tC1), -tC4*tU1, -tC4*tU2, (2.0*tC4*(tEint + tQ+1.0/tC1))/tT-2.0*tC4*tCp } };                            	 
                             
                            adAdY( 1 ) = {
                            	{                                   -tC3*tU1,                                           -tC4,          0.0,                                     (2.0*tC4*tU1)/tT },
                            	{                                 -tC3*tU1sq,                                   -2.0*tC4*tU1,          0.0,                                   (2.0*tC4*tU1sq)/tT },
                            	{                               -tC3*tU1*tU2,                                       -tC4*tU2,     -tC4*tU1,                                 (2.0*tC4*tU1*tU2)/tT },
                            	{ tC1*tCp*tU1 - tC3*tU1*(tEint + tQ+1.0/tC1), tC2*tCp - tC4*tU1sq - tC4*(tEint + tQ+1.0/tC1), -tC4*tU1*tU2, -tU1*(2.0*tC4*tCp - (2*tC4*(tEint + tQ+1.0/tC1))/tT) } };
                             
                            adAdY( 2 ) = {
                            	{                                   -tC3*tU2,          0.0,                                           -tC4,                                     (2.0*tC4*tU2)/tT },
                            	{                               -tC3*tU1*tU2,     -tC4*tU2,                                       -tC4*tU1,                                 (2.0*tC4*tU1*tU2)/tT },
                            	{                                 -tC3*tU2sq,          0.0,                                   -2.0*tC4*tU2,                                   (2.0*tC4*tU2sq)/tT },
                            	{ tC1*tCp*tU2 - tC3*tU2*(tEint + tQ+1.0/tC1), -tC4*tU1*tU2, tC2*tCp - tC4*tU2sq - tC4*(tEint + tQ+1.0/tC1), -tU2*(2.0*tC4*tCp - (2*tC4*(tEint + tQ+1.0/tC1))/tT) } };
 
                            
                            // break for TEMP case
                            break;
                        }
                    
                        default:
                        {
                            MORIS_ASSERT( false, "fn_FEM_IWG_Compressible_NS::eval_dAdY - index for Y derivative out of bounds." );
                            break;
                        }
                    }

                    // break for 2D case
                    break;
                }

                // for 3D
                case 3 :
                {
                    // get Z-velocity
                    real tU3  = tFIVelocity->val()( 3 );  
                    real tU3sq = tU3*tU3;

                    // get the 1/2 * v^2 term
                    real tQ = 0.5 * ( tU1sq * tU1sq + tU2sq * tU2sq + tU3sq * tU3sq );

                    // assemble matrices for requested state variable derivative
                    switch ( aYind )
                    {
                        // for PRESSURE
                        case 0 :
                        {
                            adAdY( 0 ) = {
                            	{ 0.0,     0.0,     0.0,     0.0,                               -tC3 },
                            	{ 0.0,     tC1,     0.0,     0.0,                           -tC3*tU1 },
                            	{ 0.0,     0.0,     tC1,     0.0,                           -tC3*tU2 },
                            	{ 0.0,     0.0,     0.0,     tC1,                           -tC3*tU3 },
                            	{ 0.0, tC1*tU1, tC1*tU2, tC1*tU3, tC1*tCp - tC3*(tEint + tQ+1.0/tC1) } };
                             
                            adAdY( 1 ) = {
                            	{ 0.0,                                tC1,         0.0,         0.0,                                 -tC3*tU1 },
                            	{ 0.0,                        2.0*tC1*tU1,         0.0,         0.0,                               -tC3*tU1sq },
                            	{ 0.0,                            tC1*tU2,     tC1*tU1,         0.0,                             -tC3*tU1*tU2 },
                            	{ 0.0,                            tC1*tU3,         0.0,     tC1*tU1,                             -tC3*tU1*tU3 },
                            	{ 0.0, tC1*tEint + tC1*tQ + tC1*tU1sq+1.0, tC1*tU1*tU2, tC1*tU1*tU3, tU1*(tC1*tCp - tC3*(tEint + tQ+1.0/tC1)) } };
                             
                            adAdY( 2 ) = {
                            	{ 0.0,         0.0,                                tC1,         0.0,                                 -tC3*tU2 },
                            	{ 0.0,     tC1*tU2,                            tC1*tU1,         0.0,                             -tC3*tU1*tU2 },
                            	{ 0.0,         0.0,                        2.0*tC1*tU2,         0.0,                               -tC3*tU2sq },
                            	{ 0.0,         0.0,                            tC1*tU3,     tC1*tU2,                             -tC3*tU2*tU3 },
                            	{ 0.0, tC1*tU1*tU2, tC1*tEint + tC1*tQ + tC1*tU2sq+1.0, tC1*tU2*tU3, tU2*(tC1*tCp - tC3*(tEint + tQ+1.0/tC1)) } };
                             
                            adAdY( 3 ) = {
                            	{ 0.0,         0.0,         0.0,                                tC1,                                 -tC3*tU3 },
                            	{ 0.0,     tC1*tU3,         0.0,                            tC1*tU1,                             -tC3*tU1*tU3 },
                            	{ 0.0,         0.0,     tC1*tU3,                            tC1*tU2,                             -tC3*tU2*tU3 },
                            	{ 0.0,         0.0,         0.0,                        2.0*tC1*tU3,                               -tC3*tU3sq },
                            	{ 0.0, tC1*tU1*tU3, tC1*tU2*tU3, tC1*tEint + tC1*tQ + tC1*tU3sq+1.0, tU3*(tC1*tCp - tC3*(tEint + tQ+1.0/tC1)) } };

                            // break for PRESSURE case
                            break;
                        }

                        // for VX
                        case 1 :
                        {
                            adAdY( 0 ) = {
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     tC1, 0.0, 0.0, 0.0,     -tC4 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{ tC1*tU1, tC2, 0.0, 0.0, -tC4*tU1 } };
                             
                            adAdY( 1 ) = {
                            	{                                tC1,         0.0,     0.0,     0.0,                                           -tC4 },
                            	{                        2.0*tC1*tU1,     2.0*tC2,     0.0,     0.0,                                    2.0*tC4*tU1 },
                            	{                            tC1*tU2,         0.0,     tC2,     0.0,                                       -tC4*tU2 },
                            	{                            tC1*tU3,         0.0,     0.0,     tC2,                                       -tC4*tU3 },
                            	{ tC1*tEint + tC1*tQ + tC1*tU1sq+1.0, 3.0*tC2*tU1, tC2*tU2, tC2*tU3, tC2*tCp - tC4*tU1sq - tC4*(tEint + tQ+1.0/tC1) } };
                             
                            adAdY( 2 ) = {
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{     tC1*tU2,     0.0,     tC2, 0.0,     -tC4*tU2 },
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{ tC1*tU1*tU2, tC2*tU2, tC2*tU1, 0.0, -tC4*tU1*tU2 } };
                             
                            adAdY( 3 ) = {
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{     tC1*tU3,     0.0, 0.0,     tC2,     -tC4*tU3 },
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{ tC1*tU1*tU3, tC2*tU3, 0.0, tC2*tU1, -tC4*tU1*tU3 } };
                            
                            // break for VX case
                            break;
                        }

                        // for VY
                        case 2 :
                        {
                            adAdY( 0 ) = {
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     tC1, 0.0, 0.0, 0.0,     -tC4 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{ tC1*tU2, 0.0, tC2, 0.0, -tC4*tU2 } };
                             
                            adAdY( 1 ) = {
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{     tC1*tU1,     tC2,     0.0, 0.0,     -tC4*tU1 },
                            	{         0.0,     0.0,     0.0, 0.0,          0.0 },
                            	{ tC1*tU1*tU2, tC2*tU2, tC2*tU1, 0.0, -tC4*tU1*tU2 } };
                             
                            adAdY( 2 ) = {
                            	{                                tC1,     0.0,         0.0,     0.0,                                           -tC4 },
                            	{                            tC1*tU1,     tC2,         0.0,     0.0,                                       -tC4*tU1 },
                            	{                        2.0*tC1*tU2,     0.0,     2.0*tC2,     0.0,                                    2.0*tC4*tU2 },
                            	{                            tC1*tU3,     0.0,         0.0,     tC2,                                       -tC4*tU3 },
                            	{ tC1*tEint + tC1*tQ + tC1*tU2sq+1.0, tC2*tU1, 3.0*tC2*tU2, tC2*tU3, tC2*tCp - tC4*tU2sq - tC4*(tEint + tQ+1.0/tC1) } };
                             
                            adAdY( 3 ) = {
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{     tC1*tU3, 0.0,     0.0,     tC2,     -tC4*tU3 },
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{ tC1*tU2*tU3, 0.0, tC2*tU3, tC2*tU2, -tC4*tU2*tU3 } };
                            
                            // break for VY case
                            break;
                        }

                        // for VZ
                        case 3 :
                        {
                            adAdY( 0 ) = {
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     0.0, 0.0, 0.0, 0.0,      0.0 },
                            	{     tC1, 0.0, 0.0, 0.0,     -tC4 },
                            	{ tC1*tU3, 0.0, 0.0, tC2, -tC4*tU3 } };
                             
                            adAdY( 1 ) = {
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{         0.0,     0.0, 0.0,     0.0,          0.0 },
                            	{     tC1*tU1,     tC2, 0.0,     0.0,     -tC4*tU1 },
                            	{ tC1*tU1*tU3, tC2*tU3, 0.0, tC2*tU1, -tC4*tU1*tU3 } };
                             
                            adAdY( 2 ) = {
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{         0.0, 0.0,     0.0,     0.0,          0.0 },
                            	{     tC1*tU2, 0.0,     tC2,     0.0,     -tC4*tU2 },
                            	{ tC1*tU2*tU3, 0.0, tC2*tU3, tC2*tU2, -tC4*tU2*tU3 } };
                             
                            adAdY( 3 ) = {
                            	{                                tC1,     0.0,     0.0,         0.0,                                           -tC4 },
                            	{                            tC1*tU1,     tC2,     0.0,         0.0,                                       -tC4*tU1 },
                            	{                            tC1*tU2,     0.0,     tC2,         0.0,                                       -tC4*tU2 },
                            	{                        2.0*tC1*tU3,     0.0,     0.0,     2.0*tC2,                                    2.0*tC4*tU3 },
                            	{ tC1*tEint + tC1*tQ + tC1*tU3sq+1.0, tC2*tU1, tC2*tU2, 3.0*tC2*tU3, tC2*tCp - tC4*tU3sq - tC4*(tEint + tQ+1.0/tC1) } };
                            
                            // break for VZ case
                            break;
                        }

                        // for TEMP
                        case 4 :
                        {
                            adAdY( 0 ) = {
                            	{                               -tC3,      0.0,      0.0,      0.0,                                  (2.0*tC4)/tT },
                            	{                           -tC3*tU1,     -tC4,      0.0,      0.0,                              (2.0*tC4*tU1)/tT },
                            	{                           -tC3*tU2,      0.0,     -tC4,      0.0,                              (2.0*tC4*tU2)/tT },
                            	{                           -tC3*tU3,      0.0,      0.0,     -tC4,                              (2.0*tC4*tU3)/tT },
                            	{ tC1*tCp - tC3*(tEint + tQ+1.0/tC1), -tC4*tU1, -tC4*tU2, -tC4*tU3, (2.0*tC4*(tEint + tQ+1.0/tC1))/tT-2.0*tC4*tCp } };
                             
                             
                            adAdY( 1 ) = {
                            	{                                   -tC3*tU1,                                           -tC4,          0.0,          0.0,                                     (2.0*tC4*tU1)/tT },
                            	{                                 -tC3*tU1sq,                                    2.0*tC4*tU1,          0.0,          0.0,                                   (2.0*tC4*tU1sq)/tT },
                            	{                               -tC3*tU1*tU2,                                       -tC4*tU2,     -tC4*tU1,          0.0,                                 (2.0*tC4*tU1*tU2)/tT },
                            	{                               -tC3*tU1*tU3,                                       -tC4*tU3,          0.0,     -tC4*tU1,                                 (2.0*tC4*tU1*tU3)/tT },
                            	{ tC1*tCp*tU1 - tC3*tU1*(tEint + tQ+1.0/tC1), tC2*tCp - tC4*tU1sq - tC4*(tEint + tQ+1.0/tC1), -tC4*tU1*tU2, -tC4*tU1*tU3, -tU1*(2.0*tC4*tCp-(2.0*tC4*(tEint + tQ+1.0/tC1))/tT) } };
                             
                             
                            adAdY( 2 ) = {
                            	{                                   -tC3*tU2,          0.0,                                           -tC4,          0.0,                                     (2.0*tC4*tU2)/tT },
                            	{                               -tC3*tU1*tU2,     -tC4*tU2,                                       -tC4*tU1,          0.0,                                 (2.0*tC4*tU1*tU2)/tT },
                            	{                                 -tC3*tU2sq,          0.0,                                    2.0*tC4*tU2,          0.0,                                   (2.0*tC4*tU2sq)/tT },
                            	{                               -tC3*tU2*tU3,          0.0,                                       -tC4*tU3,     -tC4*tU2,                                 (2.0*tC4*tU2*tU3)/tT },
                            	{ tC1*tCp*tU2 - tC3*tU2*(tEint + tQ+1.0/tC1), -tC4*tU1*tU2, tC2*tCp - tC4*tU2sq - tC4*(tEint + tQ+1.0/tC1), -tC4*tU2*tU3, -tU2*(2.0*tC4*tCp-(2.0*tC4*(tEint + tQ+1.0/tC1))/tT) } };
                             
                             
                            adAdY( 3 ) = {
                            	{                                   -tC3*tU3,          0.0,          0.0,                                           -tC4,                                     (2.0*tC4*tU3)/tT },
                            	{                               -tC3*tU1*tU3,     -tC4*tU3,          0.0,                                       -tC4*tU1,                                 (2.0*tC4*tU1*tU3)/tT },
                            	{                               -tC3*tU2*tU3,          0.0,     -tC4*tU3,                                       -tC4*tU2,                                 (2.0*tC4*tU2*tU3)/tT },
                            	{                                 -tC3*tU3sq,          0.0,          0.0,                                    2.0*tC4*tU3,                                   (2.0*tC4*tU3sq)/tT },
                            	{ tC1*tCp*tU3 - tC3*tU3*(tEint + tQ+1.0/tC1), -tC4*tU1*tU3, -tC4*tU2*tU3, tC2*tCp - tC4*tU3sq - tC4*(tEint + tQ+1.0/tC1), -tU3*(2.0*tC4*tCp-(2.0*tC4*(tEint + tQ+1.0/tC1))/tT) } };
                            
                            // break for TEMP case
                            break;
                        }
                    
                        default:
                        {
                            MORIS_ASSERT( false, "fn_FEM_IWG_Compressible_NS::eval_dAdY - index for Y derivative out of bounds." );
                            break;
                        }
                    }
                    
                    // break for 3D case
                    break;
                }
            
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dAdY() - Number of space dimensions must be 2 or 3" );
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

        void eval_dAdY_VR( 
                std::shared_ptr< Material_Model >       aMM,  
                std::shared_ptr< Constitutive_Model >   aCM,
                Field_Interpolator_Manager            * aMasterFIManager,
                const moris::Cell< MSI::Dof_Type >    & aResidualDofTypes, 
                const Matrix< DDRMat >                & aVR,
                const uint                              aI,
                Matrix< DDRMat >                      & adAdYVR )
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
            //real tC4 = tP/(tR*tT*tT);

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
                    real tQ  = ( tU1sq + tU2sq ) / 2.0;

                    // compute pre-multiplied A-matrix deriv for requested index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            adAdYVR = {
                                {                                 -tC3*tVR4,                                          0.0,                                          0.0,                                          (tC3*(2.0*tP*tVR4 - tT*tVR1))/tT },
                                {                  tC3*(tT*tVR2 - tU1*tVR4),                     -tC3*(tP*tVR4 - tT*tVR1),                                          0.0,                    -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR4 + tT*tU1*tVR1))/tT },
                                {                  tC3*(tT*tVR3 - tU2*tVR4),                                          0.0,                     -tC3*(tP*tVR4 - tT*tVR1),                    -(tC3*(tP*tT*tVR3 - 2.0*tP*tU2*tVR4 + tT*tU2*tVR1))/tT },
                                { tC3*(tT*tU1*tVR2 - tQ*tVR4 + tT*tU2*tVR3), tC3*(tP*tT*tVR2 - tP*tU1*tVR4 + tT*tU1*tVR1), tC3*(tP*tT*tVR3 - tP*tU2*tVR4 + tT*tU2*tVR1), -(tC3*(tQ*tT*tVR1 - 2.0*tP*tQ*tVR4 + tP*tT*tU1*tVR2 + tP*tT*tU2*tVR3))/tT } };
                                                     
                            break; 
                        }

                        case 1 :
                        {
                            adAdYVR = {
                                {                                                                              tC3*(tT*tVR2 - tU1*tVR4),                                                                                                                                   -tC3*(tP*tVR4 - tT*tVR1),                                                                       0.0,                                                  -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR4 + tT*tU1*tVR1))/tT },
                                {                                                                      tC3*tU1*(2.0*tT*tVR2 - tU1*tVR4),                                                                                                           2.0*tC3*(tP*tT*tVR2 - tP*tU1*tVR4 + tT*tU1*tVR1),                                                                       0.0,                                          -(tC3*tU1*(2.0*tP*tT*tVR2 - 2.0*tP*tU1*tVR4 + tT*tU1*tVR1))/tT },
                                {                                                        tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR4),                                                                                                               tC3*(tP*tT*tVR3 - tP*tU2*tVR4 + tT*tU2*tVR1),                              tC3*(tP*tT*tVR2 - tP*tU1*tVR4 + tT*tU1*tVR1),                     -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - 2.0*tP*tU1*tU2*tVR4 + tT*tU1*tU2*tVR1))/tT },
                                { tVR2 + tC3*tQ*tT*tVR2 - tC3*tQ*tU1*tVR4 + tC3*tT*tU1sq*tVR2 + tC3*tCv*tTsq*tVR2 + tC3*tT*tU1*tU2*tVR3, (tC3*((2.0*tVR1)/tC3+2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR4+2.0*tQ*tT*tVR1 - 2.0*tP*tU1sq*tVR4+2.0*tT*tU1sq*tVR1 + 6.0*tP*tT*tU1*tVR2+2.0*tP*tT*tU2*tVR3))/2.0, tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR4 + tT*tU1*tU2*tVR1), -(tC3*(tP*tQ*tT*tVR2 - 2.0*tP*tQ*tU1*tVR4 + tP*tT*tU1sq*tVR2 + tQ*tT*tU1*tVR1 + tP*tT*tU1*tU2*tVR3))/tT } };

                            break; 
                        }

                        case 2 :
                        {
                            adAdYVR = {
                                {                                                                              tC3*(tT*tVR3 - tU2*tVR4),                                                                       0.0,                                                                                                                                   -tC3*(tP*tVR4 - tT*tVR1),                                                  -(tC3*(tP*tT*tVR3 - 2.0*tP*tU2*tVR4 + tT*tU2*tVR1))/tT },
                                {                                                        tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR4),                              tC3*(tP*tT*tVR3 - tP*tU2*tVR4 + tT*tU2*tVR1),                                                                                                               tC3*(tP*tT*tVR2 - tP*tU1*tVR4 + tT*tU1*tVR1),                     -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - 2.0*tP*tU1*tU2*tVR4 + tT*tU1*tU2*tVR1))/tT },
                                {                                                                      tC3*tU2*(2.0*tT*tVR3 - tU2*tVR4),                                                                       0.0,                                                                                                           2.0*tC3*(tP*tT*tVR3 - tP*tU2*tVR4 + tT*tU2*tVR1),                                          -(tC3*tU2*(2.0*tP*tT*tVR3 - 2.0*tP*tU2*tVR4 + tT*tU2*tVR1))/tT },
                                { tVR3 + tC3*tQ*tT*tVR3 - tC3*tQ*tU2*tVR4 + tC3*tT*tU2sq*tVR3 + tC3*tCv*tTsq*tVR3 + tC3*tT*tU1*tU2*tVR2, tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR4 + tT*tU1*tU2*tVR1), (tC3*((2.0*tVR1)/tC3+2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR4+2.0*tQ*tT*tVR1 - 2.0*tP*tU2sq*tVR4+2.0*tT*tU2sq*tVR1+2.0*tP*tT*tU1*tVR2 + 6.0*tP*tT*tU2*tVR3))/2.0, -(tC3*(tP*tQ*tT*tVR3 - 2.0*tP*tQ*tU2*tVR4 + tP*tT*tU2sq*tVR3 + tQ*tT*tU2*tVR1 + tP*tT*tU1*tU2*tVR2))/tT } };
 
                            break; 
                        }
                    
                        default:
                        {
                            // error, A-matrix index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - index for A-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices

                    // break for 2D case
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

                    // compute pre-multiplied A-matrix deriv for requested index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            adAdYVR = {
                                {                                               -tC3*tVR5,                                          0.0,                                          0.0,                                          0.0,                                                           (tC3*(2.0*tP*tVR5 - tT*tVR1))/tT },
                                {                                tC3*(tT*tVR2 - tU1*tVR5),                     -tC3*(tP*tVR5 - tT*tVR1),                                          0.0,                                          0.0,                                     -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT },
                                {                                tC3*(tT*tVR3 - tU2*tVR5),                                          0.0,                     -tC3*(tP*tVR5 - tT*tVR1),                                          0.0,                                     -(tC3*(tP*tT*tVR3 - 2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT },
                                {                                tC3*(tT*tVR4 - tU3*tVR5),                                          0.0,                                          0.0,                     -tC3*(tP*tVR5 - tT*tVR1),                                     -(tC3*(tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT },
                                { tC3*(tT*tU1*tVR2 - tQ*tVR5 + tT*tU2*tVR3 + tT*tU3*tVR4), tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1), tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1), tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1), -(tC3*(tQ*tT*tVR1 - 2.0*tP*tQ*tVR5 + tP*tT*tU1*tVR2 + tP*tT*tU2*tVR3 + tP*tT*tU3*tVR4))/tT } };

                            break; 
                        }

                        case 1 :
                        {
                            adAdYVR = {
                                {                                                                                                                tC3*(tT*tVR2 - tU1*tVR5),                                                                                                                                            -tC3*(tP*tVR5 - tT*tVR1),                                                                                 0.0,                                                                       0.0,                                                                       -(tC3*(tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT },
                                {                                                                                                        tC3*tU1*(2.0*tT*tVR2 - tU1*tVR5),                                                                                                                    2.0*tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                                                                 0.0,                                                                       0.0,                                                               -(tC3*tU1*(2.0*tP*tT*tVR2 - 2.0*tP*tU1*tVR5 + tT*tU1*tVR1))/tT },
                                {                                                                                          tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR5),                                                                                                                        tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                        tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                                                       0.0,                                          -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - 2.0*tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1))/tT },
                                {                                                                                          tC3*(tT*tU1*tVR4 + tT*tU3*tVR2 - tU1*tU3*tVR5),                                                                                                                        tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                                                 0.0,                              tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                          -(tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - 2.0*tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1))/tT },
                                { (tC3*((2.0*tVR2)/tC3+2.0*tCv*tTsq*tVR2+2.0*tQ*tT*tVR2 - 2.0*tQ*tU1*tVR5+2.0*tT*tU1sq*tVR2+2.0*tT*tU1*tU2*tVR3+2.0*tT*tU1*tU3*tVR4))/2.0, (tC3*((2.0*tVR1)/tC3+2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR5+2.0*tQ*tT*tVR1 - 2.0*tP*tU1sq*tVR5+2.0*tT*tU1sq*tVR1 + 6.0*tP*tT*tU1*tVR2+2.0*tP*tT*tU2*tVR3+2.0*tP*tT*tU3*tVR4))/2.0, tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1), tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1), -(tC3*(tP*tQ*tT*tVR2 - 2.0*tP*tQ*tU1*tVR5 + tP*tT*tU1sq*tVR2 + tQ*tT*tU1*tVR1 + tP*tT*tU1*tU2*tVR3 + tP*tT*tU1*tU3*tVR4))/tT } };

                            break; 
                        }

                        case 2 :
                        {
                            adAdYVR = {
                                {                                                                                                              tC3*(tT*tVR3 - tU2*tVR5),                                                                       0.0,                                                                                                                                                -tC3*(tP*tVR5 - tT*tVR1),                                                                       0.0,                                                                       -(tC3*(tP*tT*tVR3-2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT },
                                {                                                                                        tC3*(tT*tU1*tVR3 + tT*tU2*tVR2 - tU1*tU2*tVR5),                              tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                                                                                                            tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                                                       0.0,                                          -(tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2-2.0*tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1))/tT },
                                {                                                                                                      tC3*tU2*(2.0*tT*tVR3 - tU2*tVR5),                                                                       0.0,                                                                                                                        2.0*tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                                                       0.0,                                                               -(tC3*tU2*(2.0*tP*tT*tVR3-2.0*tP*tU2*tVR5 + tT*tU2*tVR1))/tT },
                                {                                                                                        tC3*(tT*tU2*tVR4 + tT*tU3*tVR3 - tU2*tU3*tVR5),                                                                       0.0,                                                                                                                            tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                              tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                          -(tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3-2.0*tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1))/tT },
                                { (tC3*((2.0*tVR3)/tC3+2.0*tCv*tTsq*tVR3+2.0*tQ*tT*tVR3-2.0*tQ*tU2*tVR5+2.0*tT*tU2sq*tVR3+2.0*tT*tU1*tU2*tVR2+2.0*tT*tU2*tU3*tVR4))/2.0, tC3*(tP*tT*tU1*tVR3 + tP*tT*tU2*tVR2 - tP*tU1*tU2*tVR5 + tT*tU1*tU2*tVR1), (tC3*((2.0*tVR1)/tC3+2.0*tCv*tTsq*tVR1-2.0*tP*tQ*tVR5+2.0*tQ*tT*tVR1-2.0*tP*tU2sq*tVR5+2.0*tT*tU2sq*tVR1+2.0*tP*tT*tU1*tVR2+6.0*tP*tT*tU2*tVR3+2.0*tP*tT*tU3*tVR4))/2.0, tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1), -(tC3*(tP*tQ*tT*tVR3-2.0*tP*tQ*tU2*tVR5 + tP*tT*tU2sq*tVR3 + tQ*tT*tU2*tVR1 + tP*tT*tU1*tU2*tVR2 + tP*tT*tU2*tU3*tVR4))/tT } }; 
                            break; 
                        }

                        case 3 :
                        {
                            adAdYVR = {
                                {                                                                                                                tC3*(tT*tVR4 - tU3*tVR5),                                                                       0.0,                                                                       0.0,                                                                                                                                                      -tC3*(tP*tVR5 - tT*tVR1),                                                                       -(tC3*(tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT },
                                {                                                                                          tC3*(tT*tU1*tVR4 + tT*tU3*tVR2 - tU1*tU3*tVR5),                              tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                                       0.0,                                                                                                                                  tC3*(tP*tT*tVR2 - tP*tU1*tVR5 + tT*tU1*tVR1),                                          -(tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - 2.0*tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1))/tT },
                                {                                                                                          tC3*(tT*tU2*tVR4 + tT*tU3*tVR3 - tU2*tU3*tVR5),                                                                       0.0,                              tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                                                                                                  tC3*(tP*tT*tVR3 - tP*tU2*tVR5 + tT*tU2*tVR1),                                          -(tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - 2.0*tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1))/tT },
                                {                                                                                                        tC3*tU3*(2.0*tT*tVR4 - tU3*tVR5),                                                                       0.0,                                                                       0.0,                                                                                                                              2.0*tC3*(tP*tT*tVR4 - tP*tU3*tVR5 + tT*tU3*tVR1),                                                               -(tC3*tU3*(2.0*tP*tT*tVR4 - 2.0*tP*tU3*tVR5 + tT*tU3*tVR1))/tT },
                                { (tC3*((2.0*tVR4)/tC3+2.0*tCv*tTsq*tVR4+2.0*tQ*tT*tVR4 - 2.0*tQ*tU3*tVR5+2.0*tT*tU3sq*tVR4+2.0*tT*tU1*tU3*tVR2+2.0*tT*tU2*tU3*tVR3))/2.0, tC3*(tP*tT*tU1*tVR4 + tP*tT*tU3*tVR2 - tP*tU1*tU3*tVR5 + tT*tU1*tU3*tVR1), tC3*(tP*tT*tU2*tVR4 + tP*tT*tU3*tVR3 - tP*tU2*tU3*tVR5 + tT*tU2*tU3*tVR1), (tC3*((2.0*tVR1)/tC3+2.0*tCv*tTsq*tVR1 - 2.0*tP*tQ*tVR5+2.0*tQ*tT*tVR1 - 2.0*tP*tU3sq*tVR5+2.0*tT*tU3sq*tVR1+2.0*tP*tT*tU1*tVR2+2.0*tP*tT*tU2*tVR3 + 6.0*tP*tT*tU3*tVR4))/2.0, -(tC3*(tP*tQ*tT*tVR4 - 2.0*tP*tQ*tU3*tVR5 + tP*tT*tU3sq*tVR4 + tQ*tT*tU3*tVR1 + tP*tT*tU1*tU3*tVR2 + tP*tT*tU2*tU3*tVR3))/tT } };

                            break; 
                        }
                    
                        default:
                        {
                            // error, A-matrix index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dAdY - index for A-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices

                    // break for 3D case
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
        //------------------------------------------------------------------------------

        void eval_dKdY( 
                std::shared_ptr< Property >   aPropDynamicViscosity,  
                std::shared_ptr< Property >   aPropThermalConductivity,
                Field_Interpolator_Manager  * aMasterFIManager,
                const uint                    aYind,
                moris::Cell< moris::Cell< Matrix< DDRMat > > >  & adKdY )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu =  aPropDynamicViscosity->val()( 0 );
            real tLa = -2.0 * tMu / 3.0;
            real tCh = tLa + 2.0 * tMu; 

            // set size of cells
            Matrix< DDRMat > tZeroMat( tNumSpaceDims + 2, tNumSpaceDims + 2, 0.0 );
            adKdY.resize( tNumSpaceDims );
            for( uint iDim = 0; iDim < tNumSpaceDims; iDim++ ) 
            {
                adKdY( iDim ).assign( tNumSpaceDims, tZeroMat );
            }
            
            // assemble matrices based on number of spatial dimensions
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // assemble matrices for requested state variable derivative
                    switch ( aYind )
                    {
                        // for PRESSURE
                        case 0 :
                        {
                            // do nothing, derivatives are zero

                            // break for PRESSURE case
                            break;
                        }

                        // for VX
                        case 1 :
                        {
                            adKdY( 0 )( 0 )( 3, 1 ) = tCh;
                            adKdY( 0 )( 1 )( 3, 2 ) = tLa;

                            adKdY( 1 )( 0 )( 3, 2 ) = tMu;
                            adKdY( 1 )( 1 )( 3, 1 ) = tMu;
                            
                            // break for VX case
                            break;
                        }

                        // for VY
                        case 2 :
                        {
                            adKdY( 0 )( 0 )( 3, 2 ) = tMu;
                            adKdY( 0 )( 1 )( 3, 1 ) = tMu;

                            adKdY( 1 )( 0 )( 3, 1 ) = tLa;
                            adKdY( 1 )( 1 )( 3, 2 ) = tCh;
                            
                            // break for VY case
                            break;
                        }

                        // for TEMP
                        case 3 :
                        {
                            // do nothing, derivatives are zero
                            
                            // break for TEMP case
                            break;
                        }
                    
                        default:
                        {
                            MORIS_ASSERT( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY - index for Y derivative out of bounds." );
                            break;
                        }
                    }

                    // break for 2D case
                    break;
                }

                // for 3D
                case 3 :
                {
                    // assemble matrices for requested state variable derivative
                    switch ( aYind )
                    {
                        // for PRESSURE
                        case 0 :
                        {
                            // do nothing, derivatives are zero

                            // break for PRESSURE case
                            break;
                        }

                        // for VX
                        case 1 :
                        {
                            adKdY( 0 )( 0 )( 4, 1 ) = tCh;
                            adKdY( 0 )( 1 )( 4, 2 ) = tLa;
                            adKdY( 0 )( 2 )( 4, 3 ) = tLa;

                            adKdY( 1 )( 0 )( 4, 2 ) = tMu;
                            adKdY( 1 )( 1 )( 4, 1 ) = tMu;
                            //adKdY( 1 )( 2 )( 4, 3 ) = tMu;

                            adKdY( 2 )( 0 )( 4, 3 ) = tMu;
                            //adKdY( 2 )( 1 )( 4, 2 ) = tMu;
                            adKdY( 2 )( 2 )( 4, 1 ) = tMu;
                            
                            // break for VX case
                            break;
                        }

                        // for VY
                        case 2 :
                        {
                            adKdY( 0 )( 0 )( 4, 2 ) = tMu;
                            adKdY( 0 )( 1 )( 4, 1 ) = tMu;
                            //adKdY( 0 )( 2 )( 4, 3 ) = tMu;

                            adKdY( 1 )( 0 )( 4, 1 ) = tLa;
                            adKdY( 1 )( 1 )( 4, 2 ) = tCh;
                            adKdY( 1 )( 2 )( 4, 3 ) = tLa;

                            //adKdY( 2 )( 0 )( 4, 1 ) = tMu;
                            adKdY( 2 )( 1 )( 4, 3 ) = tMu;
                            adKdY( 2 )( 2 )( 4, 2 ) = tMu;
                            
                            // break for VY case
                            break;
                        }

                        // for VZ
                        case 3 :
                        {
                            adKdY( 0 )( 0 )( 4, 3 ) = tMu;
                            //adKdY( 0 )( 1 )( 4, 2 ) = tMu;
                            adKdY( 0 )( 2 )( 4, 1 ) = tMu;

                            //adKdY( 1 )( 0 )( 4, 1 ) = tMu;
                            adKdY( 1 )( 1 )( 4, 3 ) = tMu;
                            adKdY( 1 )( 2 )( 4, 2 ) = tMu;

                            adKdY( 2 )( 0 )( 4, 1 ) = tLa;
                            adKdY( 2 )( 1 )( 4, 2 ) = tLa;
                            adKdY( 2 )( 2 )( 4, 3 ) = tCh;
                            
                            // break for VZ case
                            break;
                        }

                        // for TEMP
                        case 4 :
                        {
                            // do nothing, derivatives are zero
                            
                            // break for TEMP case
                            break;
                        }
                    
                        default:
                        {
                            MORIS_ASSERT( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY - index for Y derivative out of bounds." );
                            break;
                        }
                    }
                    
                    // break for 3D case
                    break;
                }
            
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY() - Number of space dimensions must be 2 or 3" );
                    break;
                }
            }
        }

        //------------------------------------------------------------------------------

        void eval_VL_dKdY( 
                std::shared_ptr< Property >   aPropDynamicViscosity,  
                std::shared_ptr< Property >   aPropThermalConductivity,
                Field_Interpolator_Manager  * aMasterFIManager,
                const Matrix< DDRMat >      & aVL,
                const uint                    aI,
                const uint                    aJ,
                Matrix< DDRMat >            & aVLdKdY )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

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
                    real tVL4 = aVL( 3 ); 

                    // compute pre-multiplied Kij-matrix deriv for requested i-index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    aVLdKdY = {
                                        { 0.0,                  0.0,      0.0, 0.0 },
                                        { 0.0, tVL4*(tLa + 2.0*tMu),      0.0, 0.0 },
                                        { 0.0,                  0.0, tMu*tVL4, 0.0 },
                                        { 0.0,                  0.0,      0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    aVLdKdY = {
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, tMu*tVL4, 0.0 },
                                        { 0.0, tLa*tVL4,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 0
                            break; 
                        }

                        case 1 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    aVLdKdY = {
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, tLa*tVL4, 0.0 },
                                        { 0.0, tMu*tVL4,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    aVLdKdY = {
                                        { 0.0,      0.0,                  0.0, 0.0 },
                                        { 0.0, tMu*tVL4,                  0.0, 0.0 },
                                        { 0.0,      0.0, tVL4*(tLa + 2.0*tMu), 0.0 },
                                        { 0.0,      0.0,                  0.0, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 1
                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - i-index for Kij-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different Kij-matrices -> i-index

                    // break for 2D
                    break;
                }

                // for 3D
                case 3 :
                {
                    // get last value from pre-multiplication vector
                    real tVL5 = aVL( 4 ); 

                    // compute pre-multiplied Kij-matrix deriv for requested i-index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // evaluate K11
                                    aVLdKdY = { 
										{ 0.0,                  0.0,      0.0,      0.0, 0.0 },
										{ 0.0, tVL5*(tLa + 2.0*tMu),      0.0,      0.0, 0.0 },
										{ 0.0,                  0.0, tMu*tVL5,      0.0, 0.0 },
										{ 0.0,                  0.0,      0.0, tMu*tVL5, 0.0 },
										{ 0.0,                  0.0,      0.0,      0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // evaluate K12
                                    aVLdKdY = { 
										{ 0.0,      0.0,      0.0, 0, 0.0 },
										{ 0.0,      0.0, tMu*tVL5, 0, 0.0 },
										{ 0.0, tLa*tVL5,      0.0, 0, 0.0 },
										{ 0.0,      0.0,      0.0, 0, 0.0 },
										{ 0.0,      0.0,      0.0, 0, 0.0 } };
                                        

                                    break; 
                                }

                                case 2 :
                                {
                                    // evaluate K13
                                    aVLdKdY = { 
										{ 0.0,      0.0, 0.0,      0.0, 0.0 },
										{ 0.0,      0.0, 0.0, tMu*tVL5, 0.0 },
										{ 0.0,      0.0, 0.0,      0.0, 0.0 },
										{ 0.0, tLa*tVL5, 0.0,      0.0, 0.0 },
										{ 0.0,      0.0, 0.0,      0.0, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 0
                            break; 
                        }

                        case 1 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // evaluate K21
                                    aVLdKdY = { 
										{ 0.0,      0.0,      0.0, 0.0, 0.0 },
										{ 0.0,      0.0, tLa*tVL5, 0.0, 0.0 },
										{ 0.0, tMu*tVL5,      0.0, 0.0, 0.0 },
										{ 0.0,      0.0,      0.0, 0.0, 0.0 },
										{ 0.0,      0.0,      0.0, 0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // evaluate K22
                                    aVLdKdY = { 
										{ 0.0,      0.0,                  0.0,      0.0, 0.0 },
										{ 0.0, tMu*tVL5,                  0.0,      0.0, 0.0 },
										{ 0.0,      0.0, tVL5*(tLa + 2.0*tMu),      0.0, 0.0 },
										{ 0.0,      0.0,                  0.0, tMu*tVL5, 0.0 },
										{ 0.0,      0.0,                  0.0,      0.0, 0.0 } };
 
                                    break; 
                                }

                                case 2 :
                                {
                                    // evaluate K23
                                    aVLdKdY = { 
										{ 0.0, 0.0,      0.0,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0, tMu*tVL5, 0.0 },
										{ 0.0, 0.0, tLa*tVL5,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0,      0.0, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 1
                            break; 
                        }

                        case 2 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // evaluate K31
                                    aVLdKdY = { 
										{ 0.0,      0.0, 0.0,      0.0, 0.0 },
										{ 0.0,      0.0, 0.0, tLa*tVL5, 0.0 },
										{ 0.0,      0.0, 0.0,      0.0, 0.0 },
										{ 0.0, tMu*tVL5, 0.0,      0.0, 0.0 },
										{ 0.0,      0.0, 0.0,      0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // evaluate K32
                                    aVLdKdY = { 
										{ 0.0, 0.0,      0.0,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0, tLa*tVL5, 0.0 },
										{ 0.0, 0.0, tMu*tVL5,      0.0, 0.0 },
										{ 0.0, 0.0,      0.0,      0.0, 0.0 } };
 
                                    break; 
                                }

                                case 2 :
                                {
                                    // evaluate K33
                                    aVLdKdY = { 
										{ 0.0,      0.0,      0.0,                  0.0, 0.0 },
										{ 0.0, tMu*tVL5,      0.0,                  0.0, 0.0 },
										{ 0.0,      0.0, tMu*tVL5,                  0.0, 0.0 },
										{ 0.0,      0.0,      0.0, tVL5*(tLa + 2.0*tMu), 0.0 },
										{ 0.0,      0.0,      0.0,                  0.0, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 2
                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKdY - i-index for Kij-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different Kij-matrices -> i-index

                    // break for 3D case
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
                std::shared_ptr< Property >   aPropDynamicViscosity,  
                std::shared_ptr< Property >   aPropThermalConductivity,
                Field_Interpolator_Manager  * aMasterFIManager,
                const Matrix< DDRMat >      & aVR,
                const uint                    aI,
                const uint                    aJ,
                Matrix< DDRMat >            & adKdYVR )
        {
            // get the velocity FI
            Field_Interpolator * tFIVelocity =  aMasterFIManager->get_field_interpolators_for_type( { MSI::Dof_Type::VX } );
            
            // get number of Space dimensions
            uint tNumSpaceDims = tFIVelocity->get_number_of_fields();

            // get commonly used values
            //real tKa   =  aPropThermalConductivity->val()( 0 );
            real tMu   =  aPropDynamicViscosity->val()( 0 );
            real tLa   = -2.0 * tMu / 3.0;

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
                    // compute pre-multiplied Kij-matrix deriv for requested i-index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // K11
                                    adKdYVR = {
                                        { 0.0,                0.0,      0.0, 0.0 },
                                        { 0.0,                0.0,      0.0, 0.0 },
                                        { 0.0,                0.0,      0.0, 0.0 },
                                        { 0.0, tVR2*(tLa+2.0*tMu), tMu*tVR3, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // K12
                                    adKdYVR = {
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, tLa*tVR3, tMu*tVR2, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 0
                            break; 
                        }

                        case 1 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // K21
                                    adKdYVR = {
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, tMu*tVR3, tLa*tVR2, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // K22
                                    adKdYVR = {
                                        { 0.0,      0.0,                0.0, 0.0 },
                                        { 0.0,      0.0,                0.0, 0.0 },
                                        { 0.0,      0.0,                0.0, 0.0 },
                                        { 0.0, tMu*tVR2, tVR3*(tLa+2.0*tMu), 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 1
                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - i-index for Kij-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different Kij-matrices -> i-index

                    // break for 2D
                    break;
                }

                // for 3D
                case 3 :
                {


                    // compute pre-multiplied Kij-matrix deriv for requested i-index
                    switch ( aI )
                    {
                        case 0 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // K11
                                    adKdYVR = {
                                        { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                        { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, tVR2*(tLa+2.0*tMu), tMu*tVR3, tMu*tVR4, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // K12
                                    adKdYVR = {
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0, tLa*tVR3, tMu*tVR2, 0.0, 0.0 } };

                                    break; 
                                }

                                case 2 :
                                {
                                    // K13
                                    adKdYVR = {
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0, tLa*tVR4, 0.0, tMu*tVR2, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 0
                            break; 
                        }

                        case 1 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // K21
                                    adKdYVR = {
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                        { 0.0, tMu*tVR3, tLa*tVR2, 0.0, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // K22
                                    adKdYVR = {
                                        { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                        { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                        { 0.0, tMu*tVR2, tVR3*(tLa+2.0*tMu), tMu*tVR4, 0.0 } };

                                    break; 
                                }

                                case 2 :
                                {
                                    // K23
                                    adKdYVR = {
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0, tLa*tVR4, tMu*tVR3, 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 1
                            break; 
                        }

                        case 2 :
                        {
                            // compute pre-multiplied Kij-matrix deriv for requested j-index
                            switch ( aJ )
                            {
                                case 0 :
                                {
                                    // K31
                                    adKdYVR = {
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                        { 0.0, tMu*tVR4, 0.0, tLa*tVR2, 0.0 } };

                                    break; 
                                }

                                case 1 :
                                {
                                    // K32
                                    adKdYVR = {
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                        { 0.0, 0.0, tMu*tVR4, tLa*tVR3, 0.0 } };

                                    break; 
                                }

                                case 2 :
                                {
                                    // K33
                                    adKdYVR = {
                                        { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                        { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                        { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                        { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                        { 0.0, tMu*tVR2, tMu*tVR3, tVR4*(tLa+2.0*tMu), 0.0 } };

                                    break; 
                                }
                            
                                default:
                                {
                                    // error, Kiji-matrix j-index unknown
                                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - j-index for Kij-matrix out of bounds." );
                                    break;
                                }
                            } // end: switch for different Kij-matrices -> j-index
                            
                            // break for i = 2
                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKdY_VR - i-index for Kij-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different Kij-matrices -> i-index

                    // break for 3D
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

        void eval_VL_dKijidY( 
                std::shared_ptr< Property >       aPropDynamicViscosity,  
                std::shared_ptr< Property >       aPropThermalConductivity,
                Field_Interpolator_Manager      * aMasterFIManager,
                const Matrix< DDRMat >          & aVL,
                const uint                        aJ,
                moris::Cell< Matrix< DDRMat > > & aVLdKijidY )
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
            aVLdKijidY.resize( tNumSpaceDims + 1 );

            // check the pre-multiplication vector
            MORIS_ASSERT( aVL.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_VL_dKijidY - length of pre-multiplication vector incorrect." );      

            // assemble matrices based on 
            switch ( tNumSpaceDims )
            {
                // for 2D
                case 2 :
                {
                    // get last value from pre-multiplication vector
                    real tVL4 = aVL( 3 ); 

                    // compute pre-multiplied Kiji-matrix deriv for requested index
                    switch ( aJ )
                    {
                        case 0 :
                        {
                            // evaluate Ki1,i - for Y,x
                            aVLdKijidY( 1 ) = { 
                                { 0.0,                  0.0,      0.0, 0.0 },
                                { 0.0, tVL4*(tLa + 2.0*tMu),      0.0, 0.0 },
                                { 0.0,                  0.0, tMu*tVL4, 0.0 },
                                { 0.0,                  0.0,      0.0, 0.0 } };

                            // evaluate Ki1,i - for Y,y
                            aVLdKijidY( 2 ) = { 
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0, tLa*tVL4, 0.0 },
                                { 0.0, tMu*tVL4,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 } };

                            break; 
                        }

                        case 1 :
                        {
                            // evaluate Ki2,i - for Y,x
                            aVLdKijidY( 1 ) = { 
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0, tMu*tVL4, 0.0 },
                                { 0.0, tLa*tVL4,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 } };

                            // evaluate Ki2,i - for Y,x
                            aVLdKijidY( 2 ) = { 
                                { 0.0,      0.0,                  0.0, 0.0 },
                                { 0.0, tMu*tVL4,                  0.0, 0.0 },
                                { 0.0,      0.0, tVL4*(tLa + 2.0*tMu), 0.0 },
                                { 0.0,      0.0,                  0.0, 0.0 } };

                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKijidY - j-index for Kiji-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices
                    
                    // break for 2D case
                    break;
                }

                // for 3D
                case 3 :
                {
                    // get last value from pre-multiplication vector
                    real tVL5 = aVL( 4 ); 

                    // compute pre-multiplied Kiji-matrix deriv for requested index
                    switch ( aJ )
                    {
                        case 0 :
                        {
                            // evaluate Ki1,i - for Y,x
                            aVLdKijidY( 1 ) = { 
								{ 0.0,                  0.0,      0.0,      0.0, 0.0 },
								{ 0.0, tVL5*(tLa + 2.0*tMu),      0.0,      0.0, 0.0 },
								{ 0.0,                  0.0, tMu*tVL5,      0.0, 0.0 },
								{ 0.0,                  0.0,      0.0, tMu*tVL5, 0.0 },
								{ 0.0,                  0.0,      0.0,      0.0, 0.0 } };

                            // evaluate Ki1,i - for Y,y
                            aVLdKijidY( 2 ) = { 
								{ 0.0,      0.0,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0, tLa*tVL5, 0.0, 0.0 },
								{ 0.0, tMu*tVL5,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0,      0.0, 0.0, 0.0 } };
        
                            // evaluate Ki1,i - for Y,z
                            aVLdKijidY( 3 ) = { 
								{ 0.0,      0.0, 0.0,      0.0, 0.0 },
								{ 0.0,      0.0, 0.0, tLa*tVL5, 0.0 },
								{ 0.0,      0.0, 0.0,      0.0, 0.0 },
								{ 0.0, tMu*tVL5, 0.0,      0.0, 0.0 },
								{ 0.0,      0.0, 0.0,      0.0, 0.0 } };

                            break; 
                        }

                        case 1 :
                        {
                            // evaluate Ki2,i - for Y,x
                            aVLdKijidY( 1 ) = { 
								{ 0.0,      0.0,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0, tMu*tVL5, 0.0, 0.0 },
								{ 0.0, tLa*tVL5,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0,      0.0, 0.0, 0.0 },
								{ 0.0,      0.0,      0.0, 0.0, 0.0 } };

                            // evaluate Ki2,i - for Y,y
                            aVLdKijidY( 2 ) = { 
								{ 0.0,      0.0,                  0.0,      0.0, 0.0 },
								{ 0.0, tMu*tVL5,                  0.0,      0.0, 0.0 },
								{ 0.0,      0.0, tVL5*(tLa + 2.0*tMu),      0.0, 0.0 },
								{ 0.0,      0.0,                  0.0, tMu*tVL5, 0.0 },
								{ 0.0,      0.0,                  0.0,      0.0, 0.0 } };
        
                            // evaluate Ki2,i - for Y,z
                            aVLdKijidY( 3 ) = { 
								{ 0.0, 0.0,      0.0,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0, tLa*tVL5, 0.0 },
								{ 0.0, 0.0, tMu*tVL5,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0,      0.0, 0.0 } };

                            break; 
                        }

                        case 2 :
                        {
                            // evaluate Ki3,i - for Y,x
                            aVLdKijidY( 1 ) = { 
								{ 0.0,      0.0, 0.0,      0.0, 0.0 },
								{ 0.0,      0.0, 0.0, tMu*tVL5, 0.0 },
								{ 0.0,      0.0, 0.0,      0.0, 0.0 },
								{ 0.0, tLa*tVL5, 0.0,      0.0, 0.0 },
								{ 0.0,      0.0, 0.0,      0.0, 0.0 } };

                            // evaluate Ki3,i - for Y,y
                            aVLdKijidY( 2 ) = { 
								{ 0.0, 0.0,      0.0,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0, tMu*tVL5, 0.0 },
								{ 0.0, 0.0, tLa*tVL5,      0.0, 0.0 },
								{ 0.0, 0.0,      0.0,      0.0, 0.0 } };
        
                            // evaluate Ki3,i - for Y,z
                            aVLdKijidY( 3 ) = { 
								{ 0.0,      0.0,      0.0,                  0.0, 0.0 },
								{ 0.0, tMu*tVL5,      0.0,                  0.0, 0.0 },
								{ 0.0,      0.0, tMu*tVL5,                  0.0, 0.0 },
								{ 0.0,      0.0,      0.0, tVL5*(tLa + 2.0*tMu), 0.0 },
								{ 0.0,      0.0,      0.0,                  0.0, 0.0 } };

                            break; 
                        }
                    
                        default:
                        {
                            // error, A-matrix index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKijidY - j-index for Kiji-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different Kiji-matrices

                    // break for 3D case
                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_VL_dKijidY() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

        void eval_dKijidY_VR( 
                std::shared_ptr< Property >       aPropDynamicViscosity,  
                std::shared_ptr< Property >       aPropThermalConductivity,
                Field_Interpolator_Manager      * aMasterFIManager,
                const Matrix< DDRMat >          & aVR,
                const uint                        aJ,
                moris::Cell< Matrix< DDRMat > > & adKijidYVR )
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
            adKijidYVR.resize( tNumSpaceDims + 1 );

            // check the pre-multiplication vector
            MORIS_ASSERT( aVR.length() == tNumSpaceDims + 2, 
                    "fn_FEM_IWG_Compressible_NS::eval_dKijidY_VR - length of pre-multiplication vector incorrect." );      

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
                    // compute pre-multiplied Kiji-matrix deriv for requested index
                    switch ( aJ )
                    {
                        case 0 :
                        {
                            // evaluate Ki1,i - for Y,x
                            adKijidYVR( 1 ) = { 
                                { 0.0,                0.0,      0.0, 0.0 },
                                { 0.0,                0.0,      0.0, 0.0 },
                                { 0.0,                0.0,      0.0, 0.0 },
                                { 0.0, tVR2*(tLa+2.0*tMu), tMu*tVR3, 0.0 } };

                            // evaluate Ki1,i - for Y,y
                            adKijidYVR( 2 ) = { 
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, tMu*tVR3, tLa*tVR2, 0.0 } };

                            break; 
                        }

                        case 1 :
                        {
                            // evaluate Ki2,i - for Y,x
                            adKijidYVR( 1 ) = { 
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, tLa*tVR3, tMu*tVR2, 0.0 } };

                            // evaluate Ki2,i - for Y,y
                            adKijidYVR( 2 ) = { 
                                { 0.0,      0.0,                0.0, 0.0 },
                                { 0.0,      0.0,                0.0, 0.0 },
                                { 0.0,      0.0,                0.0, 0.0 },
                                { 0.0, tMu*tVR2, tVR3*(tLa+2.0*tMu), 0.0 } };

                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKijidY_VR - j-index for Kiji-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices
                    
                    // break for 2D case
                    break;
                }

                // for 3D
                case 3 :
                {
                    // compute pre-multiplied Kiji-matrix deriv for requested index
                    switch ( aJ )
                    {
                        case 0 :
                        {
                            // evaluate Ki1,i - for Y,x
                            adKijidYVR( 1 ) = { 
                                { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                { 0.0,                0.0,      0.0,      0.0, 0.0 },
                                { 0.0, tVR2*(tLa+2.0*tMu), tMu*tVR3, tMu*tVR4, 0.0 } };

                            // evaluate Ki1,i - for Y,y
                            adKijidYVR( 2 ) = { 
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0, tMu*tVR3, tLa*tVR2, 0.0, 0.0 } };
        
                            // evaluate Ki1,i - for Y,z
                            adKijidYVR( 3 ) = { 
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0, tMu*tVR4, 0.0, tLa*tVR2, 0.0 } };

                            break; 
                        }

                        case 1 :
                        {
                            // evaluate Ki2,i - for Y,x
                            adKijidYVR( 1 ) = { 
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0,      0.0,      0.0, 0.0, 0.0 },
                                { 0.0, tLa*tVR3, tMu*tVR2, 0.0, 0.0 } };

                            // evaluate Ki2,i - for Y,y
                            adKijidYVR( 2 ) = { 
                                { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                { 0.0,      0.0,                0.0,      0.0, 0.0 },
                                { 0.0, tMu*tVR2, tVR3*(tLa+2.0*tMu), tMu*tVR4, 0.0 } };
        
                            // evaluate Ki2,i - for Y,z
                            adKijidYVR( 3 ) = { 
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0, tMu*tVR4, tLa*tVR3, 0.0 } };

                            break; 
                        }

                        case 2 :
                        {
                            // evaluate Ki3,i - for Y,x
                            adKijidYVR( 1 ) = { 
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0,      0.0, 0.0,      0.0, 0.0 },
                                { 0.0, tLa*tVR4, 0.0, tMu*tVR2, 0.0 } };

                            // evaluate Ki3,i - for Y,y
                            adKijidYVR( 2 ) = { 
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0,      0.0,      0.0, 0.0 },
                                { 0.0, 0.0, tLa*tVR4, tMu*tVR3, 0.0 } };
        
                            // evaluate Ki3,i - for Y,z
                            adKijidYVR( 3 ) = { 
                                { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                { 0.0,      0.0,      0.0,                0.0, 0.0 },
                                { 0.0, tMu*tVR2, tMu*tVR3, tVR4*(tLa+2.0*tMu), 0.0 } };

                            break; 
                        }
                    
                        default:
                        {
                            // error, Kiji-matrix j-index unknown
                            MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKijidY_VR - j-index for Kiji-matrix out of bounds." );
                            break;
                        }
                    } // end: switch for different A-matrices

                    // break for 3D case
                    break;
                }                
            
                // unknown number of spatial dimensions
                default:
                {
                    MORIS_ERROR( false, "fn_FEM_IWG_Compressible_NS::eval_dKijidY_VR() - Number of space dimensions must be 2 or 3" );
                };
            }
        }

        //------------------------------------------------------------------------------

    } /* namespace fem */
} /* namespace moris */
